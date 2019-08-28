#include "sqlite3.h"
#include <unistd.h>
#include <string>
#include <iostream>
#include <vector>
#include <set>
#include "parallel_hashmap/phmap.h"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::set;
using phmap::flat_hash_map;


class SqliteHelper {

public:
    string filename;
    sqlite3 *DB;
    char *DB_ErrMsg = nullptr;
    bool file_exists;
    flat_hash_map<int, int> Qs;

    SqliteHelper(string filename, set<int> user_Qs);
    void close(){
        sqlite3_close(DB);
    }

    /*
     * if file does not exist, create seq1, seq2, minkmers, Q, value
     * if file exists, check the already written Qs
     * */


};