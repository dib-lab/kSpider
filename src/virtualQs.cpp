#include "virtualQs.hpp"
#include <iostream>
#include <fstream>
#include "combinations.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <sqlite3.h>
using boost::adaptors::transformed;
using boost::algorithm::join;
using std::cerr;
using std::endl;
using std::to_string;

inline uint64_t create_mask(unsigned kSize, unsigned Q) {
    return ((1ULL << Q * 2ULL) - 1ULL) << (kSize * 2ULL - Q * 2ULL);
}

virtualQs::virtualQs(string index_prefix, set<int> allQs) {

    // Load the sqlite DB
    // TODO assert file exist before trying to load.

    string sqlite_file = index_prefix + "_kCluster.sqlite";
    int rc;
    rc = sqlite3_open(sqlite_file.c_str(), &this->DB);
    if(rc) {
        fprintf(stderr, "Can't open ths sqlite database: %s\n", sqlite3_errmsg(this->DB));
        exit(0);
    }

    // Loading the index
    cerr << "[INFO] Loading the index" << endl;
    this->KF = kDataFrame::load(index_prefix);
    this->kSize = KF->ksize();
    this->index_prefix = index_prefix;

    // Constructing masks
    for (auto const & Q : allQs)
        this->mainQs.insert(Q);

    bool ksize_in_Qs = (this->masks.find(this->kSize) != this->masks.end());

    if (!ksize_in_Qs)
        this->mainQs.insert(this->kSize);

    for (auto const &Q : this->mainQs) {
        this->masks[Q] = create_mask(this->kSize, Q);
        this->superColors[Q] = flat_hash_map<uint64_t, flat_hash_set<uint64_t>>();
        this->superColorsCount[Q] = flat_hash_map<uint64_t, uint64_t>();
        this->temp_superColors[Q] = flat_hash_set<uint64_t>();
    }

    string colors_map = this->index_prefix + "colors.intvectors";
    ifstream input(colors_map.c_str());
    int size;
    input >> size;
    color_to_ids = flat_hash_map<uint64_t, std::vector<uint32_t>>(size);
    for (int i = 0; i < size; i++) {
        uint64_t color, colorSize;
        input >> color >> colorSize;
        uint32_t sampleID;
        color_to_ids[color] = std::vector<uint32_t>(colorSize);
        for (int j = 0; j < colorSize; j++) {
            input >> sampleID;
            color_to_ids[color][j] = sampleID;
        }
    }


    // Read NamesMap

//    ifstream namesMapIn(index_prefix + ".namesMap");
//    namesMapIn >> size;
//    for (int i = 0; i < size; i++) {
//        uint32_t sample_id;
//        string sample_name;
//        namesMapIn >> sample_id >> sample_name;
//        this->namesMap[sample_id] = sample_id;
//    }

}

uint64_t virtualQs::create_super_color(flat_hash_set<uint64_t> &colors) {
    uint64_t seed = colors.size();
    for (auto &color : colors)
        seed ^= color + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

void virtualQs::calculate_kmers_number() {

    for (auto const &superColor : this->superColors[this->kSize]) {
        vector<uint32_t> tr_ids;
        for (auto const &color : superColor.second)
            for (auto const &id : this->color_to_ids[color])
                tr_ids.push_back(id);

        uint32_t color_count = this->superColorsCount[this->kSize][superColor.first];

        for (auto const &tr_id : tr_ids) {
            bool tr_exist = (this->seq_to_kmers_no.find(tr_id) != this->seq_to_kmers_no.end());
            if (tr_exist) {
                this->seq_to_kmers_no[tr_id] += color_count;
            } else {
                this->seq_to_kmers_no[tr_id] = color_count;
            }
        }

    }

}

void virtualQs::pairwise() {
    this->calculate_kmers_number();
    Combo combo = Combo();
    cerr << "Processing Qs: ";
    for (auto const &Q : this->mainQs) {
        cerr << Q << " ";
        for (auto const &superColor : this->superColors[Q]) {
            vector<uint32_t> tr_ids;
            for (auto const &color : superColor.second)
                for (auto const &id : this->color_to_ids[color])
                    tr_ids.push_back(id);

            uint32_t color_count = this->superColorsCount[Q][superColor.first];
            combo.combinations(tr_ids.size());
            for (auto const &seq_pair : combo.combs) {
                uint32_t _seq1 = tr_ids[seq_pair.first];
                uint32_t _seq2 = tr_ids[seq_pair.second];
                _seq1 > _seq2 ? this->edges[{{_seq1, _seq2}, Q}] += color_count : this->edges[{{_seq2, _seq1}, Q}] += color_count;
            }
        }

    }
    cerr << endl;
}


inline string prepare_insertion(string & Q_names, vector<uint32_t> &values) {
    std::string VALUES = join(values | transformed([](uint32_t d) { return std::to_string(d); }), ",");
    return "INSERT INTO virtualQs (seq1, seq2, min_kmers, " + Q_names + ") VALUES (" + VALUES + ");";
}


void virtualQs::export_to_tsv(){
    Combo combo = Combo();
    combo.combinations(this->seq_to_kmers_no.size());
    std::ofstream myfile;
    myfile.open(this->index_prefix + "_kCluster.tsv");
    myfile << "ID" << '\t' << "seq1" << '\t' << "seq2" << '\t' << "min_kmers" << '\t';
    string delimiter = "";
    for(auto const &val : this->mainQs){
        myfile << delimiter << "Q_" << val;
        delimiter = '\t';
    }
    myfile << '\n';

    uint64_t line_count = 0;


    for (auto const &seq_pair : combo.combs) {
        uint32_t seq1;
        uint32_t seq2;
        if(seq_pair.first > seq_pair.second){
            seq1 = seq_pair.first + 1;
            seq2 = seq_pair.second + 1;
        }else{
            seq2 = seq_pair.first + 1;
            seq1 = seq_pair.second + 1;
        }

        uint32_t min_kmers = std::min(this->seq_to_kmers_no[seq1], this->seq_to_kmers_no[seq2]);

        vector<uint32_t> values;
        uint32_t sum = 0;
        for(auto const & Q : this->mainQs){
            uint32_t val = this->edges[{{seq1, seq2},Q}];
            sum += val;
            values.emplace_back(val);
        }

        if(sum){
            myfile << line_count << '\t' << seq1 << '\t' << seq2 << '\t' << min_kmers << '\t';
            delimiter = "";
            for(auto const &val : values){
                myfile << delimiter << val;
                delimiter = '\t';
            }
            myfile << '\n';
        }
        line_count ++;
    }
    myfile.close();
}

void virtualQs::export_to_sqlite() {
    this->superColors.clear();
    this->superColorsCount.clear();
    vector<uint32_t > values;
    int rc;
    values.reserve(4);
    sqlite3_exec(this->DB, "PRAGMA cache_size=10000000", NULL, NULL, &this->DB_ErrMsg);
    sqlite3_exec(this->DB, "BEGIN TRANSACTION", NULL, NULL, &this->DB_ErrMsg);
    sqlite3_exec(this->DB, "PRAGMA synchronize = OFF", NULL, NULL, &this->DB_ErrMsg);
    sqlite3_exec(this->DB, "PRAGMA jorunal_mode = MEMORY", NULL, NULL, &this->DB_ErrMsg);
    std::string Q_names = "Q_" + join(this->mainQs | transformed([](uint8_t d) { return std::to_string(d); }), ", Q_");

    Combo combo = Combo();
    combo.combinations(this->seq_to_kmers_no.size());

    int total_elements = 3 + this->mainQs.size();

    for (auto const &seq_pair : combo.combs) {
        uint32_t seq1;
        uint32_t seq2;
        if(seq_pair.first > seq_pair.second){
            seq1 = seq_pair.first + 1;
            seq2 = seq_pair.second + 1;
        }else{
            seq2 = seq_pair.first + 1;
            seq1 = seq_pair.second + 1;
        }

        uint32_t min_kmers = std::min(this->seq_to_kmers_no[seq1], this->seq_to_kmers_no[seq2]);

        vector<uint32_t> values = {seq1, seq2, min_kmers};
        uint32_t sum = 0;
        for(auto const & Q : this->mainQs){
            uint32_t val = this->edges[{{seq1, seq2},Q}];
            sum += val;
            values.emplace_back(val);
        }

        if(sum){
//            cout << "inserting: " << prepare_insertion(Q_names, values) << endl;
            rc = sqlite3_exec(this->DB, prepare_insertion(Q_names, values).c_str(), NULL, 0, &this->DB_ErrMsg);
            if( rc != SQLITE_OK ){
                fprintf(stderr, "SQL error: %s\n", this->DB_ErrMsg);
                sqlite3_free(&this->DB_ErrMsg);
            }
        }

    }


    sqlite3_exec(this->DB, "END TRANSACTION", NULL, NULL, &this->DB_ErrMsg);
    sqlite3_close(this->DB);
}