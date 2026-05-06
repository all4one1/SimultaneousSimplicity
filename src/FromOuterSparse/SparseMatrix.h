#pragma once
#include <vector>
#include <map>
#include <string>
// #define DEBUG

// format of storage:
// Compressed Sparse Row
struct SparseMatrix
{
	enum side
	{
		center,
		east,
		north,
		west,
		south,
		front,
		back
	};

	int Nfull = 0;	// the input (linear) size of a matrix
	int nval = 0;	// number of non-zero elements
	int nrow = 0;	// number of raws
	double zero_threshold = 1e-15;

	std::vector <double> val;
	std::vector <int> col;
	std::vector <int> row;
	std::vector <double> diag;
	std::vector <int> type;


	SparseMatrix();
	SparseMatrix(int n_);
	~SparseMatrix();


	void add_one_next(double v, int i, int t = -1);
	void insert_one(int place_in_vector, int matrix_row, double value, int column, int t = -1);
	void erase_one(int place_in_vector, int matrix_row);
	void add_line_with_map(std::map<int, double> elements, int current_line);
	void endline(int l);


	void make_sparse_2d_laplace(int nx, int ny, double a = 1);
	void make_sparse_from_double_array(int n_, double** M);
	void make_sparse_from_2d_vector(std::vector<std::vector<double>> M);
	double get_element(int ii, int jj);
	int get_index(int ii, int jj);
	int get_type(int ii, int jj);
	void set_type(int ii, int jj, int t);

	void update_diag();
	double max_element_abs();

	double get_diag(int l);
	
	void update(int ii, int jj, double value);
	void resize(int n_);

	void reset();

	SparseMatrix& operator=(std::vector<std::vector<double>> SM)
	{
		this->make_sparse_from_2d_vector(SM);
		return *this;
	};
	double& operator()(int ii, int jj);

	// implementaion of [][] operator
	struct Bracket
	{
		SparseMatrix& root;
		int i;
		Bracket(SparseMatrix& r, int i_) : root(r), i(i_) {}
		double& operator[](int j) { return root(i, j); }
	};
	Bracket operator[] (int i)
	{
		return Bracket(*this, i);
	}


	//TODO: erase zeros after some manipulations (like a decomposition technique)
	void erase_zeros();

	// return S = sum_j (Aij*yj), i.e. sum up all elements in a line (row)
	double line(int q, double* y);

	// auxiliary sum for LU solver
	double line1(int q, double* y);
	// auxiliary sum for LU solver
	double line2(int q, double* y);


	// IO defined in the second *.cpp file:
	void save_compressed_matrix(std::string filename = "CompressedMatrix.dat");
	void save_compressed_matrix_with_rhs(double* b, std::string filename = "CompressedMatrix.dat");
	void save_compressed_matrix_binary(std::string filename = "CompressedMatrix.bin");
	void read_compressed_matrix(std::string filename = "CompressedMatrix.dat");
	void read_compressed_matrix_with_rhs(double* b, std::string filename = "CompressedMatrix.dat");
	void read_compressed_matrix_binary(std::string filename = "CompressedMatrix.bin");

	void save_full_matrix(int precision = 4, std::string filename = "FullMatrix.dat");
	void save_full_matrix_v2(int precision = 4, std::string filename = "FullMatrix.dat");
	void save_full_matrix_with_rhs(int precision, double* b, std::string filename = "FullMatrix.dat");
	void read_full_matrix_with_rhs(int N, double* b, std::string filename = "FullMatrix.dat");

	void print_storage();
	void save_type();
	void print_index_ij(int l);
	void print_compressed_matrix();
	void print_full_matrix();

};
