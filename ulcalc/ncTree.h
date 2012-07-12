/*
 *  ncTree.h
 *
 *  Created by Christoph Nägeli <christoph.naegeli@psi.ch> on 02.07.12.
 *
 */

#ifndef NCTREE_H
#define NCTREE_H

#include <iostream>
#include <vector>

template<class T>
class ncTree {
	
	public:
		ncTree(T node);
		~ncTree();
		
		ncTree *getChild1();
		ncTree *getChild2();
		T* getData() {return &fData;}
		
		void dump(bool parantheses = true);
		std::string toString(bool parantheses = true);
		
		void setChild1(ncTree *newC);
		void setChild2(ncTree *newC);
	private:
		ncTree<T> *child1;
		ncTree<T> *child2;
		T fData; // data associated to the children
};

template<class T>
ncTree<T>::ncTree(T node) : child1(NULL), child2(NULL) {
	fData = node;
} // ncTree()

template<class T>
ncTree<T>::~ncTree()
{
	if (child1) delete child1;
	if (child2) delete child2;
} // ~ncTree()

template<class T>
ncTree<T>* ncTree<T>::getChild1() {return child1;}

template<class T>
ncTree<T>* ncTree<T>::getChild2() {return child2;}

template<class T>
void ncTree<T>::setChild1(ncTree<T> *newC) {child1 = newC;}

template<class T>
void ncTree<T>::setChild2(ncTree<T> *newC) {child2 = newC;}

template<class T>
void ncTree<T>::dump(bool parantheses)
{
	std::cout << toString(parantheses) << std::endl;
} // dump()

template<class T>
std::string ncTree<T>::toString(bool parantheses)
{
	std::string result;
	
	if (child1) {
		if(parantheses) result.append("(");
		result.append(child1->toString(parantheses));
		if(parantheses) result.append(")");
	}
	result.append(fData.toString());
	if (child2) {
		if(parantheses) result.append("(");
		result.append(child2->toString(parantheses));
		if(parantheses) result.append(")");
	}
	
	return result;
} // toString()

#endif
