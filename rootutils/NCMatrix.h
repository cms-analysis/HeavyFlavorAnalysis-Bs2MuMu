/*
 *  NCMatrix.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 08.11.12.
 *
 */

#include <vector>

#include <TGraph.h>

class NCMatrix {
	
	typedef std::vector<double> row_t;
	typedef std::string out_of_bounds_t;
	
	public:
		NCMatrix() {}
		
		double getElement(int row, int col);
		void setElement(int row, int col, double value);
		
		size_t getNRows();
		size_t getNCols();
		
		void print();
		void clear();
		
		void readFile(const char *name);
		
		// Visualisation
		TGraph *plot(int xrow, int yrow);
	
	public:
		static NCMatrix *matrixFromFile(const char *name);
	
	private:
		bool checkRange(int row, int col);
	
	private:
		std::vector<row_t> fMatrix;
};
