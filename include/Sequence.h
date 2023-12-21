#ifndef Sequence_h_
#define Sequence_h_

#include <string>
#include <vector>
using namespace std;

class Sequence {
public:
    vector<string> seqVec;
    string geneName;
    long numSeq;
    Sequence(string _geneName): geneName(_geneName), numSeq(0) {}
    Sequence(vector<string> _seqVec, string _name);
    Sequence(): geneName("NA"), numSeq(0) {}

    void setGeneName(string _geneName) {
        geneName = _geneName;
    }

    string getGeneName() {
        return geneName;
    }

    string getSeq(long i) {
        return seqVec[i];
    }

    long getNumSeq() {
        return numSeq;
    }

    void addSeq(string given_seq);
    void checkAddSeq(string given_seq, string name);
};

#endif