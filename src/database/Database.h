#pragma once
#include "sqlite3.h"
#include <string>
#include <vector>
#include <map>
#include <iostream> 
#include <sstream> 
#include <iomanip>
using std::string;
using std::vector; 
using std::cout;
using std::endl;
using std::to_string;




struct Database
{
    sqlite3* db = nullptr, *back = nullptr;
    sqlite3_backup* backup = nullptr;
    char* zErrMsg = 0;
    int rc, read_n = 0;
    string query;

    sqlite3_stmt* stmt_read_blob = nullptr, *stmt_write_blob = nullptr;
    sqlite3_stmt* stmt_read_double = nullptr, * stmt_write_double = nullptr;
    sqlite3_stmt* st1 = NULL, * st2 = NULL, * st3 = NULL;
    sqlite3_stmt* stmt_read = NULL, * stmt_write = NULL, *stmt_line = NULL;
    sqlite3_stmt* stmt_insert = NULL;
    Database(string name);
    Database() {};
    void open_database(string name);
    void close_database();

    void create_table(string table, vector<string> one_val, string primary);
    void replace_by_str(std::string table, std::string pKey, std::string pKeyVal, std::string colName, double val);

    void make_query(string q);
    void insert(string table, string key, string val);
    void insert_many(string table, string primary_key, string primary_key_value, std::map<string, double> col_val);
    void insert(string table, std::pair<string, string> p);
    void init_rows(string table, string key, std::vector<std::string> val, std::map<string, double> m = {});
    void update(string table, string key_name, string key_value, string column, string value);

    void update(string table, string key, string key_value, std::map<string, double> m);


    void write_blob(string table, string field_name, string key, double *data, unsigned int Nbytes);
    void read_blob(string table, string field_name, string key, double** data, unsigned int Nbytes = 0);
   
    int read_int(string table, string column, string key);
    double read_double(string table, string primary_key, string key, string column);

    double read_double_ps(string table, string param);
    void write_double_ps(string name, double val);

    void test_back(const std::string& source_file, const std::string& destination_file)
    {
         // Открываем исходную базу данных
        //int rc = sqlite3_open(source_file.c_str(), &source_db);
        //if (rc != SQLITE_OK) {
        //    std::cerr << "Failed to open source database: " << sqlite3_errmsg(source_db) << std::endl;
        //    sqlite3_close(source_db);
        //    return;
        //}

        // Открываем целевую базу данных
        //rc = sqlite3_open(destination_file.c_str(), &destination_db);
        //if (rc != SQLITE_OK) {
        //    std::cerr << "Failed to open destination database: " << sqlite3_errmsg(destination_db) << std::endl;
        //    sqlite3_close(source_db);
        //    sqlite3_close(destination_db);
        //    return;
        //}

        //rc = sqlite3_open("test.db", &db);

        rc = sqlite3_open("back.db", &back);

       // source_db = db;
       // destination_db = back;

        // Создаем резервную копию
        backup = sqlite3_backup_init(back, "main", db, "main");
        if (backup) {
            // Копируем данные
            rc = sqlite3_backup_step(backup, -1);  // Копируем все страницы
            if (rc != SQLITE_DONE) {
                std::cerr << "Failed to copy database: " << sqlite3_errmsg(back) << std::endl;
            }
            else {
                std::cout << "Database backup completed successfully!" << std::endl;
            }
            sqlite3_backup_finish(backup);  // Завершаем копирование
        }
        else {
            std::cerr << "Failed to create backup: " << sqlite3_errmsg(back) << std::endl;
        }
    }
};

