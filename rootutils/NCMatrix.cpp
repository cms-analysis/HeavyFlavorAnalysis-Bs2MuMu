/*
 *  NCMatrix.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 08.11.12.
 *
 */

#include "NCMatrix.h"

#include <iostream>

#include <TString.h>

NCMatrix *NCMatrix::matrixFromFile(const char *name)
{
	NCMatrix *m = new NCMatrix;
	m->readFile(name);
	
	return m;
} // matrixFromFile

bool NCMatrix::checkRange(int row, int col)
{
	return (0 <= row && (size_t)row < getNRows()) && (0 <= col && (size_t)col < getNCols());
} // checkRange()

double NCMatrix::getElement(int row, int col)
{
	double value = 0.0;
	
	if (!checkRange(row, col))
		throw out_of_bounds_t(Form("NCMatrix: row=%d,col=%d exceeding matrix",row,col));
	
	if ((size_t)col < fMatrix[row].size())
		value = fMatrix[row][col];
	
	return value;
} // getElement()

void NCMatrix::setElement(int row, int col, double value)
{
	while (row >= (int)fMatrix.size())
		fMatrix.push_back(row_t());
	
	while (col >= (int)fMatrix[row].size())
		fMatrix[row].push_back(0.0);
	
	fMatrix[row][col] = value;
} // setElement()

size_t NCMatrix::getNRows()
{
	return fMatrix.size();
} // getNRows()

size_t NCMatrix::getNCols()
{
	size_t ncol = 0;
	for (size_t j = 0; j < fMatrix.size(); j++) {
		if (fMatrix[j].size() > ncol)
			ncol = fMatrix[j].size();
	}
	
	return ncol;
} // getNCols()

void NCMatrix::print()
{
	using std::cout; using std::endl;
	size_t ncols = getNCols(), nrows = getNRows();
	size_t r,c;
	
	for(r = 0; r < nrows; r++) {
		for (c = 0; c < ncols; c++)
			cout << '\t' << getElement(r, c);
		cout << endl;
	}
} // print()

void NCMatrix::clear()
{
	fMatrix.clear();
} // clear()

void NCMatrix::readFile(const char *name)
{
	FILE *file = fopen(name,"r");
	char *buffer;
	size_t len;
	std::string line;
	std::string delim(" \t,\n");
	std::string::iterator b,e;
	double value;
	row_t *row;
	
	clear();
	
	if (!file) {
		std::cerr << "Unable to open file '" << name << "'" << std::endl;
		goto bail;
	}
	
	while( (buffer = fgetln(file, &len)) != NULL) {
		row = NULL;
		if (buffer[0] == '#')
			continue; // comment line
		
		line = std::string(buffer,len);
		
		b = line.begin();
		do {
			e = find_first_of(b, line.end(), delim.begin(), delim.end());
			
			if ((e - b) > 0) {
				if (!row) {
					fMatrix.push_back(row_t());
					row = &(fMatrix.back());
				}
				value = atof((std::string(line,b - line.begin(), (e-b))).c_str());
				row->push_back(value);
			}
			
			b = e + 1;
		} while (e != line.end());
	}
	
bail:
	fclose(file);
} // readFile()

TGraph *NCMatrix::plot(int xcol, int ycol)
{
	TGraph *result = new TGraph;
	int nrows = getNRows(),j;
	
	checkRange(0, xcol);
	checkRange(0, ycol);
	
	result->SetName(Form("x%dy%d",xcol,ycol));
	result->SetTitle(Form("x-axis col %d, y-axis col %d",xcol,ycol));
	
	for(j = 0; j < nrows; j++)
		result->SetPoint(j, getElement(j, xcol), getElement(j, ycol));
	
	return result;
} // plotRow()
