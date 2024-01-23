#include "func/Sequence.h"

Sequence::Sequence(vector<string> _seqVec, string _name) {
    seqVec = _seqVec;
    geneName = _name;
    numSeq = seqVec.size();
}

void Sequence::addSeq(string given_seq) {
    seqVec.push_back(given_seq);
    numSeq += 1;
}


void Sequence::checkAddSeq(string given_seq, string name) {
    if (geneName == name) {
        addSeq(given_seq);
    }
}