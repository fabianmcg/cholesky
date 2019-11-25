#ifndef __ELIMITATION_TREE_MATRIX_OPERATIONS_H__
#define __ELIMITATION_TREE_MATRIX_OPERATIONS_H__

#include <vector>

#include "../matrix-formats/matrix-formats.h"
#include "../data-structures/data-structures.h"

namespace __core__ {
namespace __linear_algebra__ {
namespace __matrix_operations__ {
template <typename T,typename IT,typename Allocator,SparseCXSFormat frmt> LRTree<IT,IT,Allocator,int> elimitationTree(const MatrixCXS<T,IT,Allocator,frmt>& matrix) {
	LRTree<IT,IT,Allocator,int> tree(matrix.n()+64,-1);
	std::vector<int> vars(matrix.n(),-1);
	Node<IT,IT,int> *node,*nodeRow;
	for(IT row=1;row<matrix.n();++row) {
		for(IT k=matrix.ptr(row);row<matrix.ptr(row+1);++k) {
			IT col=matrix.indxs(k);
			if(col>=row)
				break;
			if(vars[row]==-1) {
				nodeRow=&(tree.insert(row,tree.invalidPos,tree.nilNode));
				vars[row]=nodeRow->self;
			}
			nodeRow=&(tree[vars[row]]);
			if(vars[col]==-1) {
				node=&(tree.insert(col,row,*nodeRow));
				vars[col]=node->self;
			}
			else {
				node=&(tree[vars[col]]);
				IT ancestor=node->value;
				node->value=row;
				bool update=(ancestor!=row);
				while((ancestor!=tree.invalidPos)&&update) {
					node=&(tree[vars[ancestor]]);
					ancestor=node->value;
					node->value=row;
					if(ancestor==row) {
						update=false;
						break;
					}
				}
				if(update)
					tree.moveUnder(node->self,nodeRow->self);
			}
		}
	}
	for(size_t row=0;row<matrix.n();++row)
		if(vars[row]==-1)
			tree.insert(row,tree.invalidPos,tree.nilNode);
	return tree;
}

template <typename T,typename IT,typename Allocator> std::vector<IT> descendentCount(LRTree<IT,IT,Allocator,int>& tree) {

}
template <typename T,typename IT,typename Allocator> std::vector<IT> postOrderTree(LRTree<IT,IT,Allocator,int>& tree) {
	std::vector<IT> pi(tree.size());
	IT count=0,tmp=0,parent=0,off=0;
	auto reorder=[&pi,&count,&tmp,&parent,&off](Node<IT,IT,int>& node,long depth) {
		node.key=((parent+tmp)-count)+off;

		++count;
	};
}
}
}
}
#endif
