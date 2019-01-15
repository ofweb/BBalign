package main

import (
	"bytes"
	"fmt"
	"log"
	"strings"
	"sync"
	"time"

	"go.uber.org/atomic"
	"go.uber.org/zap"
	"go.uber.org/zap/zapcore"
)

type template struct {
	sequence []rune //actual dna sequence
	//TODO add
	//  UPAC ambiguty codes
	//  planned dna mutations
	//  support for CODING reqions (CDS)
	//  codontables and translation
	//  planned aa mutations
}

func (t *template) MarshalLogObject(oe zapcore.ObjectEncoder) error {
	oe.AddString("sequence", string(t.sequence))
	return nil
}

type seqRes struct {
	sequence []rune //actual dna sequence
	//TODO: add
	//  quality array
	//  multiple seqres inputs
}

func (s *seqRes) MarshalLogObject(oe zapcore.ObjectEncoder) error {
	oe.AddString("sequence", string(s.sequence))
	return nil
}

//we store where we skip ahead on either of the sequences
//this lets us store a lot less per alignemnt
// but getting the actual alignment is a bit more complicated
type offset struct {
	offset     int
	onTemplate bool
	offsetLoc  int
	otherLoc   int
}

func (o offset) MarshalLogObject(oe zapcore.ObjectEncoder) error {
	oe.AddInt("offset", o.offset)
	oe.AddBool("ontemplate", o.onTemplate)
	oe.AddInt("offsetLoc", o.offsetLoc)
	oe.AddInt("otherLoc", o.otherLoc)
	return nil
}

type offsets []offset

func (os offsets) MarshalLogArray(oa zapcore.ArrayEncoder) error {
	for _, o := range os {
		oa.AppendObject(o)
	}
	return nil
}

//TODO a way to get some nice mutations
// add classic ascii visual
// add json output

//edits made, we can change this to be an ummtable list github.com/mndrix/ps as each child node only addes offsets

func (os offsets) pushOffset(log *zap.SugaredLogger, o offset) offsets {
	if len(os) == 0 {
		return append(os, o)
	}

	for _, of := range os {
		if of.offsetLoc+of.offset > o.offsetLoc {
			log.Panicw("trying to push an offset that is before the last in the list",
				"list", os,
				"new offset", o,
				"conflict", of,
			)
		}
	}

	last := &os[len(os)-1]
	if last.onTemplate == o.onTemplate {
		if o.offsetLoc == last.offsetLoc+last.offset {
			log.Infow("mergeing offset",
				"list", os,
				"offset", o,
			)
			last.offset = last.offset + o.offset
			return os
		}
	}
	log.Infow("appending offset",
		"list", os,
		"offset", o,
	)
	return append(os, o)
}

type alignment struct {
	offsets     offsets
	score       int
	templatePos int //last position in template sequence
	seqResPos   int //last position in seqRes sequence
}

func (a alignment) MarshalLogObject(oe zapcore.ObjectEncoder) error {
	oe.AddArray("offsets", a.offsets)
	oe.AddInt("score", a.score)
	oe.AddInt("templatePos", a.templatePos)
	oe.AddInt("seqResPos", a.seqResPos)
	return nil
}

func (a *alignment) asciiVisual(log *zap.SugaredLogger, t *template, sr *seqRes) string {
	var templateOut bytes.Buffer
	var seqResOut bytes.Buffer
	var opLine bytes.Buffer
	for ti, si := 0, 0; ti < len(t.sequence) && si < len(sr.sequence) && ti <= a.templatePos && si <= a.seqResPos; {
		in, offset := inOffset(a, ti, si)
		if in {
			if offset.onTemplate {
				templateOut.WriteString(string(t.sequence[ti : ti+offset.offset]))
				seqResOut.WriteString(strings.Repeat("-", offset.offset))
				opLine.WriteString(strings.Repeat("-", offset.offset))
				ti += offset.offset
			} else {
				seqResOut.WriteString(string(sr.sequence[si : si+offset.offset]))
				templateOut.WriteString(strings.Repeat("-", offset.offset))
				opLine.WriteString(strings.Repeat("-", offset.offset))
				si += offset.offset
			}
			continue
		}
		if t.sequence[ti] == sr.sequence[si] {
			seqResOut.WriteString(string(sr.sequence[si]))
			templateOut.WriteRune('|')
			opLine.WriteString(string(t.sequence[ti]))
		} else {
			seqResOut.WriteString(string(sr.sequence[si]))
			templateOut.WriteRune('.')
			opLine.WriteString(string(t.sequence[ti]))
		}
		ti++
		si++
	}
	templateOut.WriteRune('\n')
	templateOut.Write(opLine.Bytes())
	templateOut.WriteRune('\n')
	templateOut.Write(seqResOut.Bytes())
	return templateOut.String()

}

type scoring interface {
	upperBound(*alignment, *template, *seqRes) int
	lowerBound(*alignment, *template, *seqRes) int
	score(*alignment, *template, *seqRes) int
}

type bestAlignment struct {
	mutex *sync.RWMutex
	a     *alignment
}

type node struct {
	running  *atomic.Int32 //current nodes where execution have started
	spawned  *atomic.Int32 //total nodes spawned
	finished *atomic.Int32 //total nodes finished
	active   *atomic.Int32 //current nodes that are not finished
	t        *template
	s        *seqRes
	pa       alignment //partial alignmnet
	sf       scoring
	cba      *bestAlignment //current best alignment
}

func (n node) MarshalLogObject(oe zapcore.ObjectEncoder) error {
	oe.AddObject("partial Alignment", n.pa)
	return nil
}

func (n *node) isComplete() bool {
	return n.isDoneWithSeqRes() && n.isdoneWithTemplate()
}

func (n *node) isdoneWithTemplate() bool {
	return n.pa.templatePos >= (len(n.t.sequence) - 1)
}

func (n *node) isDoneWithSeqRes() bool {
	return n.pa.seqResPos >= (len(n.s.sequence) - 1)
}

func (n *node) updateCBA(log *zap.SugaredLogger) {
	score := n.sf.score(&n.pa, n.t, n.s)
	n.pa.score = score
	n.cba.mutex.RLock()
	bscore := n.cba.a.score
	n.cba.mutex.RUnlock()
	log.Infow("end node",
		"node", n,
		"current best", bscore,
	)
	if score < bscore {
		log.Infow("new best alignment",
			"score", score,
			"current best", bscore,
		)
		//our alignment is better
		n.cba.mutex.Lock()
		//check that cba did not change while we were waiting for the lock
		if score < n.cba.a.score {
			n.cba.a = &n.pa
		}
		n.cba.mutex.Unlock()
	}
}

func processNode(log *zap.SugaredLogger, n node) {
	n.running.Add(1)
	defer n.running.Dec()
	defer n.active.Dec()
	defer n.finished.Add(1)

	log.Infow("staring on node",
		"node", n,
	)

	//alignment is complete so report that
	if n.isComplete() {
		n.updateCBA(log)
		return
	}

	//we are out of sequence on the sequence result, so skip the rest of the template and report
	if n.isDoneWithSeqRes() {
		offset := offset{
			offset:     len(n.t.sequence) - 1 - n.pa.templatePos,
			onTemplate: true,
			offsetLoc:  n.pa.templatePos,
		}
		log.Infow("done with seqres",
			"node", n,
			"new offset", offset,
		)
		n.pa.offsets = n.pa.offsets.pushOffset(log, offset)
		n.pa.templatePos = len(n.t.sequence) - 1
		n.updateCBA(log)
		return
	}

	//we are out of sequence on the template, so skip the rest of the sequencing result and report
	if n.isdoneWithTemplate() {
		offset := offset{
			offset:     len(n.s.sequence) - 1 - n.pa.seqResPos,
			onTemplate: false,
			offsetLoc:  n.pa.seqResPos,
		}
		log.Infow("done with seqres",
			"node", n,
			"new offset", offset,
		)
		n.pa.offsets = n.pa.offsets.pushOffset(log, offset)
		n.pa.seqResPos = len(n.s.sequence) - 1
		n.updateCBA(log)
		return
	}

	n.cba.mutex.RLock()
	bscore := n.cba.a.score
	n.cba.mutex.RUnlock()

	// we handle three cases, match/mismatch ins or del
	matchNode := n
	matchNode.pa.templatePos++
	matchNode.pa.seqResPos++
	log.Infow("new match node",
		"node", matchNode,
	)
	if bscore > matchNode.sf.lowerBound(&matchNode.pa, matchNode.t, matchNode.s) {
		go processNode(log, matchNode)
		n.spawned.Add(1)
		n.active.Add(1)
		//log.Infow("adding matchnode")
	}

	insNode := n
	insOffset := offset{
		offset:     1,
		onTemplate: true,
		offsetLoc:  insNode.pa.templatePos,
	}
	insNode.pa.offsets = insNode.pa.offsets.pushOffset(log, insOffset)
	insNode.pa.templatePos++
	log.Infow("insnode",
		"node", insNode,
	)
	if bscore > insNode.sf.lowerBound(&insNode.pa, insNode.t, insNode.s) {
		go processNode(log, insNode)
		n.spawned.Add(1)
		n.active.Add(1)
		//log.Infow("adding ins node")
	}

	delNode := n
	delOffset := offset{
		offset:     1,
		onTemplate: false,
		offsetLoc:  delNode.pa.seqResPos,
	}
	delNode.pa.offsets = delNode.pa.offsets.pushOffset(log, delOffset)
	delNode.pa.seqResPos++
	log.Infow("delnode",
		"node", delNode,
	)
	if bscore > delNode.sf.lowerBound(&delNode.pa, delNode.t, delNode.s) {
		go processNode(log, delNode)
		n.spawned.Add(1)
		n.active.Add(1)
		//log.Infow("adding del node")
	}

	return
}

type simpleScoring struct {
	match    int
	mismatch int
	skip     int
}

func (s *simpleScoring) upperBound(a *alignment, t *template, sr *seqRes) int {
	//all mismatch + skip for all overflow
	score := s.score(a, t, sr)
	tempMissing := len(t.sequence) - a.templatePos
	seqResMissing := len(sr.sequence) - a.seqResPos
	if tempMissing >= seqResMissing {
		score += (tempMissing - seqResMissing) * s.skip
		score += seqResMissing * s.mismatch
	} else {
		score += (seqResMissing - tempMissing) * s.skip
		score += tempMissing * s.mismatch
	}
	return score
}
func (s *simpleScoring) lowerBound(a *alignment, t *template, sr *seqRes) int {
	//all match + skip for all overflow
	score := s.score(a, t, sr)
	tempMissing := len(t.sequence) - a.templatePos
	seqResMissing := len(sr.sequence) - a.seqResPos
	if tempMissing >= seqResMissing {
		score += (tempMissing - seqResMissing) * s.skip
		score += seqResMissing * s.match
	} else {
		score += (seqResMissing - tempMissing) * s.skip
		score += tempMissing * s.match
	}
	return score
}
func (s *simpleScoring) score(a *alignment, t *template, sr *seqRes) int {
	//only for completed part for alignment
	score := 0
	ti := 0
	si := 0
	if a.seqResPos == 0 {
		return s.skip * len(sr.sequence)
	}
	if a.templatePos == 0 {
		return s.skip * len(t.sequence)
	}

	for ti < len(t.sequence) && si < len(sr.sequence) && ti <= a.templatePos && si <= a.seqResPos {
		in, offset := inOffset(a, ti, si)
		if in {
			score += s.skip * offset.offset
			if offset.onTemplate {
				ti += offset.offset
			} else {
				si += offset.offset
			}
			continue
		}
		if t.sequence[ti] == sr.sequence[si] {
			score += s.match
		} else {
			score += s.mismatch
		}
		ti++
		si++
	}
	return score
}

func inOffset(a *alignment, ti int, si int) (bool, offset) {
	for _, offset := range a.offsets {
		if offset.onTemplate && (ti >= offset.offsetLoc && ti <= offset.offsetLoc+offset.offset) {
			return true, offset
		}
		if !offset.onTemplate && (si >= offset.offsetLoc && si <= offset.offsetLoc+offset.offset) {
			return true, offset
		}
	}
	return false, offset{}
}

func main() {
	logger, err := zap.NewDevelopment()
	if err != nil {
		log.Fatalf("can't initialize zap logger: %v", err)
	}
	defer logger.Sync()
	sLogger := logger.Sugar()

	temp := template{sequence: []rune("TAGTA")}
	contig := seqRes{sequence: []rune("AGTA")}

	a := alignment{
		offsets:     make([]offset, 0),
		score:       65535,
		templatePos: 0,
		seqResPos:   0,
	}

	ba := bestAlignment{
		mutex: &sync.RWMutex{},
		a:     &a,
	}

	scoring := simpleScoring{match: 1, mismatch: 5, skip: 5}

	running := atomic.NewInt32(0)
	spawned := atomic.NewInt32(0)
	finished := atomic.NewInt32(0)
	active := atomic.NewInt32(1) //1 for the root that we will spawn

	node := node{
		t:        &temp,
		s:        &contig,
		sf:       &scoring,
		pa:       a,
		cba:      &ba,
		running:  running,
		spawned:  spawned,
		finished: finished,
		active:   active,
	}

	go processNode(sLogger, node)

	start := time.Now()
	//prevActive := active.Load()
	for {
		a := active.Load()

		if a == 0 {
			//output the alignment
			break
		}

		ba.mutex.RLock()
		bscore := ba.a.score
		ba.mutex.RUnlock()

		/* 		if a == prevActive {
			sLogger.Fatalw("no change",
				"running", running.Load(),
				"spawned", spawned.Load(),
				"finished", finished.Load(),
				"active", active.Load(),
				"runtime", time.Since(start),
				"current best", bscore,
			)
		} */

		//prevActive = a

		sLogger.Infow("not done yet",
			"running", running.Load(),
			"spawned", spawned.Load(),
			"finished", finished.Load(),
			"active", active.Load(),
			"runtime", time.Since(start),
			"current best", bscore,
		)

		if time.Since(start).Seconds() > float64(20) {
			sLogger.Fatalw("out of time",
				"running", running.Load(),
				"spawned", spawned.Load(),
				"finished", finished.Load(),
				"active", active.Load(),
				"runtime", time.Since(start),
				"current best", bscore,
			)
		}

		if a >= 1000 {
			sLogger.Fatalw("too many gone-rouge-teens",
				"running", running.Load(),
				"spawned", spawned.Load(),
				"finished", finished.Load(),
				"active", active.Load(),
				"runtime", time.Since(start),
				"current best", bscore,
			)
		}

		time.Sleep(time.Microsecond * 2)
	}

	//we are done
	sLogger.Sync()
	ba.mutex.RLock()
	fmt.Printf("alignment with score: %d\n", ba.a.score)
	fmt.Printf("%+v\n", ba.a)
	fmt.Println(ba.a.asciiVisual(sLogger, &temp, &contig))
}
