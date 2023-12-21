#include "util/util.h"


long encode_nt_to_num(string nt) {
    long ans = 0;
    for (size_t i = 0; i < nt.size(); i++) {
        char a = nt[i];
        long num;
        switch (a) {
            case 'A':
                num = 0;
                break;
            case 'C':
                num = 1;
                break;
            case 'G':
                num = 2;
                break;
            case 'T':
                num = 3;
                break;
            default:
                num = -1;
                break;
        }
        ans = ans * 4 + num;
    }
    return ans;
}
