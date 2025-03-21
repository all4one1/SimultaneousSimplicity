// TODO: skip if all zeros

#include "Database.h"
#include <iostream> 
#include <sstream> 
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::to_string;

Database::Database(string name)
{
    open_database(name);  
    //open_database("back", back);
}
void Database::open_database(string name)
{
    name.append(".db");
    rc = sqlite3_open(name.c_str(), &db);
    
    if (rc)
    {
        cout << stderr << " Can't open database: " << sqlite3_errmsg(db) << endl;
    }
    else {
        //cout << "Opened database successfully" << endl;
    }
};
void Database::close_database()
{
    sqlite3_close(db);
}

std::vector<string> lastCol, lastArg;
static int callback(void* NotUsed, int argc, char** argv, char** azColName) {
    int i;

    if (argc > 0)
    {
        lastCol.clear();
        lastArg.clear();
    }

    for (i = 0; i < argc; i++) 
    {
         lastCol.push_back(azColName[i]);
         lastArg.push_back(argv[i] ? argv[i] : "NULL");
    }
    return 0;
}

void Database::create_table(string table, std::vector<string> one_val, string primary)
{
    query = "CREATE TABLE IF not exists " + table + "(";

    for (auto& it : one_val)
    {
        query.append(it + ",");
    }
    query.append("PRIMARY KEY(" + primary + "))");
    make_query(query);
}
void Database::make_query(string q)
{
    //q = "select count(*) from pragma_table_info('fields')";
    rc = sqlite3_exec(db, q.c_str(), callback, 0, &zErrMsg);
   
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    else {
       // fprintf(stdout, "Query made successfully\n");
    }
}


void Database::insert(string table, string key, string val)
{
    string query = "INSERT OR IGNORE INTO " + table + "(" + key + ") VALUES(" + val + ")";
    make_query(query);
}
void Database::insert_many(string table, string primary_key, string primary_key_value, std::map<string, double> col_val)
{
    auto str = [](double d, int prec = 10)
    {
        std::ostringstream s;
        s << std::setprecision(prec) << d;
        return s.str();
    };
    insert(table, primary_key, primary_key_value);

    query = "UPDATE " + table + " SET ";
    for (auto it = col_val.begin(); it != col_val.end(); ++it)
    {
        query += it->first + " = " + str(it->second);
        if (it != std::prev(col_val.end()))
            query += ", ";
    }
    query += " WHERE " + primary_key + " = " + primary_key_value;
    make_query(query);
}
void Database::insert(string table, std::pair<string, string> p)
{
    string query = "INSERT OR IGNORE INTO " + table + "(" + p.first + ") VALUES(" + p.second + ")";
    make_query(query);
}
void Database::init_rows(string table, string key, std::vector<string> val, std::map<string, double> m)
{
    for (auto& it : val) 
    {
        insert(table, key, it);
        if (!m.empty())
            update(table, "field", key, m);
    }

}
void Database::update(string table, string key_name, string key_value, string column, string value)
{
    string query = "UPDATE " + table + " SET " + column + " = " + value + " WHERE " + key_name + " = " + key_value;
    rc = sqlite3_exec(db, query.c_str(), callback, 0, &zErrMsg);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    else {
       // fprintf(stdout, "Records created successfully\n");
    }
}
void Database::update(string table, string key, string key_value, std::map<string, double> m)
{
    auto str = [](double d, int prec = 10)
    {
        std::ostringstream s;
        s << std::setprecision(prec) << d;
        return s.str();
    };
    
    string query = "UPDATE " + table + " SET ";
    
    for (auto it = m.begin(); it != m.end(); ++it)
    {
       // query += it->first + " = " + std::to_string(it->second);
        query += it->first + " = " + str(it->second);
        if (it != std::prev(m.end()))
            query += ", ";
    }
    query += " WHERE " + key + " = " + key_value;

    make_query(query);
}

void Database::replace_by_str(std::string table, std::string pKey, std::string pKeyVal, std::string colName, double val)
{
    query = "REPLACE INTO " + table + "(" + pKey + "," + colName + ") VALUES ('" + pKeyVal + "'," + to_string(val) + ")";
    make_query(query);
}
void Database::write_blob(string table, string field_name, string key, double *data, unsigned int Nbytes)
{
    //update or insert?
    if (stmt_write_blob == nullptr)
    {
        string query = "UPDATE " + table + " SET data = :data WHERE field = :field";
        //string query = "REPLACE INTO " + table + "(field, data) VALUES(:field, :data)";
        rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt_write_blob, NULL);
        if (rc != SQLITE_OK)    cout << "prepare failed: " << sqlite3_errmsg(db) << endl;
    }
    rc = sqlite3_bind_blob(stmt_write_blob, sqlite3_bind_parameter_index(stmt_write_blob, ":data"), data, Nbytes, SQLITE_STATIC);
    rc = sqlite3_bind_text(stmt_write_blob, sqlite3_bind_parameter_index(stmt_write_blob, ":field"), key.c_str(), -1, SQLITE_STATIC);

    if (rc != SQLITE_OK)
    {
        cout << "bind failed: " << sqlite3_errmsg(db) << endl;
    }
    else
    {
        rc = sqlite3_step(stmt_write_blob);
        if (rc != SQLITE_DONE)  cout << "execution failed: " << sqlite3_errmsg(db) << endl;
    }
    rc = sqlite3_reset(stmt_write_blob);
}
void Database::read_blob(string table, string field_name, string key, double **data, unsigned int Nbytes)
{
    if (stmt_read_blob == nullptr)
    {
        string query = "SELECT data from " + table + " WHERE " + field_name + " = ? ";
        rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt_read_blob, NULL);
        if (rc != SQLITE_OK)    cout << "prepare failed: " << sqlite3_errmsg(db) << endl;
    }
   // rc = sqlite3_bind_blob(st2, 1, data, Nbytes, SQLITE_STATIC);
    rc = sqlite3_bind_text(stmt_read_blob, 1, key.c_str(), -1, SQLITE_STATIC);
    rc = sqlite3_step(stmt_read_blob);

    if (rc == SQLITE_ROW)
    {
        Nbytes = sqlite3_column_bytes(stmt_read_blob, 0);
        memcpy(*data, sqlite3_column_blob(stmt_read_blob, 0), Nbytes);
    }
    else
    {
        std::cout << "execution failed: " << sqlite3_errmsg(db) << endl;
    }
    rc = sqlite3_reset(stmt_read_blob);
}


int Database::read_int(string table, string column, string key)
{
    query = "SELECT " + column + " from parameters WHERE id = " + key;
    rc = sqlite3_exec(db, query.c_str(), callback, 0, &zErrMsg);

    int val = -1;
    if (column == lastCol.back())
    {
        std::istringstream(lastArg.back()) >> val;
    }

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    else {
     //   fprintf(stdout, "Query made successfully\n");
    }

    return val;
}
double Database::read_double(string table, string primary_key, string key, string column)
{
    query = "SELECT " + column + " from parameters WHERE " + primary_key + " = " + key;
    cout << query << endl;
    rc = sqlite3_exec(db, query.c_str(), callback, 0, &zErrMsg);

    double val = NAN;
    if (column == lastCol.back())
    {
        std::istringstream(lastArg.back()) >> val;
    }

    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    else {
      //  fprintf(stdout, "Query made successfully\n");
    }

    return val;
}

double Database::read_double_ps(string table, string param)
{
    double res = NAN;
    if (stmt_read_double == nullptr)
    {
        string query = "SELECT value from " + table + " WHERE name = ? ";
        rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt_read_double, NULL);
        if (rc != SQLITE_OK)    cout << "prepare failed: " << sqlite3_errmsg(db) << endl;
    }
    // rc = sqlite3_bind_blob(st2, 1, data, Nbytes, SQLITE_STATIC);
    rc = sqlite3_bind_text(stmt_read_double, 1, param.c_str(), -1, SQLITE_STATIC);
    rc = sqlite3_step(stmt_read_double);

    if (rc == SQLITE_ROW)
    {
        res = sqlite3_column_double(stmt_read_double, 0);
    }
    else
    {
        std::cout << "execution failed: " << sqlite3_errmsg(db) << endl;
    }
    rc = sqlite3_reset(stmt_read_double);
    return res;
}
void Database::write_double_ps(string name, double val)
{
    if (stmt_write_double == NULL)
    {
        string query = "UPDATE Parameters SET value = :value WHERE name = :name";
        //string query = "REPLACE INTO " + table + "(field, data) VALUES(:field, :data)";
        rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt_write_double, NULL);
        if (rc != SQLITE_OK)    cout << "prepare failed: " << sqlite3_errmsg(db) << endl;
    }

    rc = sqlite3_bind_double(stmt_write_double, 1, val);
    rc = sqlite3_bind_text(stmt_write_double, 2, name.c_str(), -1, SQLITE_STATIC);


    if (rc != SQLITE_OK)
    {
        cout << "bind failed: " << sqlite3_errmsg(db) << endl;
    }
    else
    {
        rc = sqlite3_step(stmt_write_double);
        if (rc != SQLITE_DONE)  cout << "execution failed: " << sqlite3_errmsg(db) << endl;
    }
    rc = sqlite3_reset(stmt_write_double);
}
