#include "sqlite_helper.hpp"


SqliteHelper::SqliteHelper(string filename, set<int> user_Qs) {
    this->filename = filename;

    int rc;

    for (auto const &Q : user_Qs)
        this->Qs[Q]++;

    file_exists = access(filename.c_str(), F_OK) != -1;

    if (!file_exists) { // sqlite file does not exist
        cerr << "creating the database" << endl;
        rc = sqlite3_open(filename.c_str(), &this->DB);
        if (rc) {
            fprintf(stderr, "Can't open ths sqlite database: %s\n", sqlite3_errmsg(this->DB));
            exit(0);
        }

        string sql_statement = "CREATE TABLE virtualQs" \
                               "(ID INTEGER PRIMARY KEY AUTOINCREMENT," \
                               "seq1            INT     NOT NULL," \
                               "seq2            INT     NOT NULL," \
                               "min_kmers       INT     NOT NULL," \
                               "Q       INT     NOT NULL," \
                               "shared       INT     NOT NULL" \
                               ");";

        rc = sqlite3_exec(this->DB, sql_statement.c_str(), NULL, 0, &DB_ErrMsg);
        if( rc != SQLITE_OK ){
            fprintf(stderr, "SQL error: %s\n", DB_ErrMsg);
            sqlite3_free(&DB_ErrMsg);
        }
    }
}


