/*
 *  ncFormula.cpp
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 02.07.12.
 *
 */

// my headers
#include "ncFormula.h"

// Standard headers
#include <algorithm>
#include <utility>
#include <queue>

using std::string;

ncFormula::ncFormula() : fName(""), fOp(false)
{}

ncFormula::ncFormula(std::string name, bool op) : fName(name), fOp(op)
{}

#pragma mark -

static ncTree<ncFormula> *read_exp1(string::iterator b, string::iterator e,string::iterator &last);

// read single expression
static ncTree<ncFormula> *read_exp6(string::iterator b, string::iterator e, string::iterator &last)
{
	ncTree<ncFormula> *result;
	
	if (*b == '(') {
		int paraCount = 1;
		last = b+1;
		while (paraCount > 0 && last != e) {
			if (*last == '(')		paraCount++;
			else if (*last == ')')	paraCount--;
			last++;
		}
		if (paraCount > 0) throw string("Paranthesis does not close");
		result = read_exp1(b+1, last-1, last);
		last++; // skipp paranthesis
	} else if (*b == '!') {
		result = new ncTree<ncFormula>(ncFormula("!",true));
		result->setChild2(read_exp6(b+1,e,last));
	} else {
		char specials[] = {'_','.'}; // append here?
		for (last = b; last != e &&
			 ( ('a' <= *last && *last <= 'z') ||
			  ('A' <= *last && *last <= 'Z') ||
			  ('0' <= *last && *last <= '9') ||
			  (std::find(specials, specials+sizeof(specials)/sizeof(char), *last) != specials+sizeof(specials)/sizeof(char)) ); ++last) {}
		result = new ncTree<ncFormula>(ncFormula(string(b,last),false));
	}
	
	return result;
} // read_exp6()

// parse *,/
static ncTree<ncFormula> *read_exp5(string::iterator b, string::iterator e, string::iterator &last)
{
	ncTree<ncFormula> *result = read_exp6(b, e, last);
	
	while (last != e && (*last == '*' || *last == '/') ) {
		
		ncTree<ncFormula> *child1 = result;
		result = new ncTree<ncFormula>(ncFormula(string(last,last+1),true));
		b = last+1;
		ncTree<ncFormula> *child2 = read_exp6(b, e, last);
		
		result->setChild1(child1);
		result->setChild2(child2);
	}
	
	return result;
} // read_exp5()

// parse +,-
static ncTree<ncFormula> *read_exp4(string::iterator b, string::iterator e, string::iterator &last)
{
	ncTree<ncFormula> *result = read_exp5(b,e,last);
	
	while (last != e && (*last == '+' || *last == '-') ) {
		ncTree<ncFormula> *child1 = result;
		result = new ncTree<ncFormula>(ncFormula(string(last,last+1),true));
		b = last+1;
		ncTree<ncFormula> *child2 = read_exp5(b, e, last);
		
		result->setChild1(child1);
		result->setChild2(child2);
	}
	
	return result;
} // read_exp4()

// parse <=,...
static ncTree<ncFormula> *read_exp3(string::iterator b, string::iterator e, string::iterator &last)
{
	ncTree<ncFormula> *result = read_exp4(b,e,last);
	char ops [] = {'<','>','=','!'};
	char *ops_end = ops + sizeof(ops)/sizeof(char);
	string opname;
	
	// skip the operator
	do {
		opname = "";
		while (last != e && (std::find(ops, ops_end, *last) != ops_end))
			opname.push_back(*last++);
		
		if (opname.size() > 0) {
			b = last;
			ncTree<ncFormula> *child1 = result;
			ncTree<ncFormula> *child2 = read_exp4(b,e,last);
			
			result = new ncTree<ncFormula>(ncFormula(opname,true));
			result->setChild1(child1);
			result->setChild2(child2);
		}
	} while (opname.size() > 0);
	
	return result;
} // read_exp3()

// parse &&
static ncTree<ncFormula> *read_exp2(string::iterator b, string::iterator e, string::iterator &last)
{
	string opName = "&&";
	ncTree<ncFormula> *result = read_exp3(b,e,last);
	
	while (last != e && std::search(last,e,opName.begin(),opName.end()) == last) {
		b = last + 2;
		ncTree<ncFormula> *child1 = result;
		ncTree<ncFormula> *child2 = read_exp3(b, e, last);
		
		result = new ncTree<ncFormula>(ncFormula(opName,true));
		result->setChild1(child1);
		result->setChild2(child2);
	}
	
	return result;
} // read_exp2()

// parse ||
static ncTree<ncFormula> *read_exp1(string::iterator b, string::iterator e, string::iterator &last)
{
	string opName = "||";
	ncTree<ncFormula> *result = read_exp2(b, e, last);
	
	while (last != e && std::search(last,e,opName.begin(),opName.end()) == last) {
		b = last + 2;
		ncTree<ncFormula> *child1 = result;
		ncTree<ncFormula> *child2 = read_exp2(b, e, last);
		
		result = new ncTree<ncFormula>(ncFormula(opName,true));
		result->setChild1(child1);
		result->setChild2(child2);
	}
	
	return result;
} // read_exp1()

#pragma mark -

ncTree<ncFormula> *read_formula(std::string formula)
{
	std::string::iterator it;
	
	it = std::remove(formula.begin(),formula.end(),' ');
	it = std::remove(formula.begin(),it,'\t');
	it = std::remove(formula.begin(),it,'\n');
	it = std::remove(formula.begin(),it,'\r');
	
	formula.erase(it,formula.end());
	
	return read_exp1(formula.begin(),formula.end(),it);
} // read_formula()

std::set<string> get_cuts(ncTree<ncFormula> *tree)
{
	std::set<string> result;
	std::queue<ncTree<ncFormula> *> q;
	ncFormula *formula;
	
	// iterate the tree
	q.push(tree);
	while (!q.empty()) {
		tree = q.front();
		formula = tree->getData();
		if ( formula->isOp() && (formula->getName().compare("||") == 0 || formula->getName().compare("&&") == 0) ) {
			if (tree->getChild1())
				q.push(tree->getChild1());
			if (tree->getChild2())
				q.push(tree->getChild2());
		} else
			result.insert(tree->toString(false));
		
		q.pop();
	}
	
	return result;
} // get_cuts()

std::set<string> get_dependencies(ncTree<ncFormula> *tree)
{
	std::set<string> result;
	std::queue<ncTree<ncFormula> *> q;
	ncFormula *formula;
	
	// iterate the tree
	q.push(tree);
	while (!q.empty()) {
		tree = q.front();
		formula = tree->getData();
		if (formula->isOp()) {
			if (tree->getChild1()) q.push(tree->getChild1());
			if (tree->getChild2()) q.push(tree->getChild2());
		} else {
			string s = formula->getName();
			if(std::find_if(s.begin(),s.end(),isalpha) != s.end())
				result.insert(s);
		}
		q.pop();
	}
	
	return result;
} // get_dependencies()
