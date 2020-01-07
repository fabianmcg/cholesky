#ifndef __SUPER_NODES_H__
#define __SUPER_NODES_H__

#include "../data-structures/data-structures.h"

namespace __core__ {
namespace __linear_algebra__ {
namespace __cholesky__ {
struct SuperNode {
	int parent=-1;
	int evars=0;
	std::vector<int> nodes=std::vector<int>();
	SuperNode(int v) {
	}
	SuperNode(int p,int v,int ev=0) :nodes(std::vector<int>()) {
		parent=p;
		nodes.push_back(v);
		evars=ev;
	}
	void insert(int v) {
		nodes.push_back(v);
	}
};
namespace __super_nodes__ {
template <typename IT,typename Allocator> void simplifyTree(LRTree<IT,IT,Allocator,int>& tree,size_t branchSize,double tolerance=0.2) {
	typedef Node<IT,IT,int> NT;
	IT bbsize=branchSize+tolerance*branchSize;
	bool down=true;
	int node=(*tree).left_child;
	size_t i=0;
	while((i++)<tree.size()) {
		if(down) {
			if(tree[node].value<bbsize) {
				tree.eraseChildren(tree[node],tree[node].value);
				down=false;
			}
			else {
				size_t evars=1;
				while(evars<branchSize) {
					IT it=tree[node].right_child;
					while((it!=tree.invalidPos)&&(evars<branchSize)) {
						NT n=tree[it];
						tree.erase(n);
						it=n.left_sibling;
						++evars;
					}
				}
				down=true;
			}
		}
		if(down&&(tree[node].left_child!=tree.invalidPos))
			node=tree[node].left_child;
		else if(tree[node].right_sibling!=tree.invalidPos) {
			node=tree[node].right_sibling;
			down=true;
		}
		else {
			node=tree[node].parent;
			down=false;
		}
		if(node==0||node==tree.invalidPos)
			break;
	}
}
template <typename IT,typename Allocator> void populateTree(LRTree<IT,SuperNode,Allocator,int>& sntree,LRTree<IT,IT,Allocator,int>& simplifiedTree) {
	size_t pos=1;
	int node=(*simplifiedTree).left_child;
	bool down=true;
	size_t i=0;
	Node<IT,SuperNode,int>* n=&(*sntree);
	while((i++)<simplifiedTree.size()) {
		if(down) {
			n=&sntree.insert(pos++,SuperNode(n->key,simplifiedTree[node].self,simplifiedTree[node].value+1),*n);
			while(simplifiedTree[node].left_child!=simplifiedTree.invalidPos) {
				node=simplifiedTree[node].left_child;
				n=&sntree.insert(pos++,SuperNode(n->key,simplifiedTree[node].self,simplifiedTree[node].value+1),*n);
			}
		}
		down=false;
		if(simplifiedTree[node].right_sibling!=simplifiedTree.invalidPos) {
			node=simplifiedTree[node].right_sibling;
			n=&(sntree[n->parent]);
			down=true;
		}
		else {
			node=simplifiedTree[node].parent;
			n=&(sntree[n->parent]);
		}
		if(node==0||node==simplifiedTree.invalidPos)
			break;
	}
}
}
std::ostream& operator<<(std::ostream& ost,const SuperNode& node){
	ost<<node.evars<<" "<<node.parent<<" {";
	for(size_t i=0;i<node.nodes.size();++i) {
		if(i!=(node.nodes.size()-1))
			ost<<node.nodes[i]<<", ";
		else
			ost<<node.nodes[i];
	}
	ost<<"}";
	return ost;
}
template <typename IT,typename Allocator> auto superNodes(LRTree<IT,IT,Allocator,int>& tree,size_t branchSize,double tolerance=0.2) {
	LRTree<IT,IT,Allocator,int> simplifiedTree(tree);
	__super_nodes__::simplifyTree(simplifiedTree,branchSize,tolerance);
	LRTree<IT,SuperNode,Allocator,int> sntree(simplifiedTree.size(),-1);
	__super_nodes__::populateTree(sntree,simplifiedTree);
//	simplifiedTree.print(std::cerr)<<std::endl;
//	sntree.print(std::cerr)<<std::endl;
	int node=(*sntree).left_child;
	bool down=true;
	size_t i=0;
	while((i++)<sntree.size()) {
		if(down)
			while(sntree[node].left_child!=sntree.invalidPos)
				node=sntree[node].left_child;
		down=false;
		if((sntree[node].left_sibling==sntree.invalidPos)&&(sntree[node].right_sibling==sntree.invalidPos)&&(sntree[node].left_child==sntree.invalidPos)&&(sntree[node].parent!=0)) {
			Node<IT,SuperNode,int>& tmp=sntree[node];
			node=tmp.parent;
			sntree.erase(tmp);
		}
		else {
			if(sntree[node].right_sibling!=sntree.invalidPos) {
				size_t evars=sntree[node].value.evars;
				int it=sntree[node].right_sibling;
				while(it!=sntree.invalidPos) {
					if((evars+sntree[it].value.evars)<=branchSize) {
						sntree[node].value.evars+=sntree[it].value.evars;
						sntree[node].value.nodes.insert(sntree[node].value.nodes.end(),sntree[it].value.nodes.begin(),sntree[it].value.nodes.end());
						evars+=sntree[it].value.evars;
						int tmp=sntree[it].right_sibling;
						sntree.erase(sntree[it]);
						it=tmp;
						continue;
					}
					it=sntree[it].right_sibling;
				}
				if(sntree[node].right_sibling!=sntree.invalidPos) {
					node=sntree[node].right_sibling;
					down=true;
				}
			}
			else
				node=sntree[node].parent;
		}
		if(node==0||node==sntree.invalidPos)
			break;
	}
	return sntree;
}
template <typename IT,typename Allocator> auto superNodes(LRTree<IT,IT,Allocator,int>& tree,LRTree<IT,IT,Allocator,int>& simplifiedTree,size_t branchSize,double tolerance=0.2) {
	simplifiedTree.import(tree);
	__super_nodes__::simplifyTree(simplifiedTree,branchSize,tolerance);
	LRTree<IT,SuperNode,Allocator,int> sntree(simplifiedTree.size(),-1);
	__super_nodes__::populateTree(sntree,simplifiedTree);
//	simplifiedTree.print(std::cerr)<<std::endl;
//	sntree.print(std::cerr)<<std::endl;
	int node=(*sntree).left_child;
	bool down=true;
	size_t i=0;
	while((i++)<sntree.size()) {
		if(down)
			while(sntree[node].left_child!=sntree.invalidPos)
				node=sntree[node].left_child;
		down=false;
		if((sntree[node].left_sibling==sntree.invalidPos)&&(sntree[node].right_sibling==sntree.invalidPos)&&(sntree[node].left_child==sntree.invalidPos)&&(sntree[node].parent!=0)) {
			Node<IT,SuperNode,int>& tmp=sntree[node];
			node=tmp.parent;
			sntree.erase(tmp);
		}
		else {
			if(sntree[node].right_sibling!=sntree.invalidPos) {
				size_t evars=sntree[node].value.evars;
				int it=sntree[node].right_sibling;
				while(it!=sntree.invalidPos) {
					if((evars+sntree[it].value.evars)<=branchSize) {
						sntree[node].value.evars+=sntree[it].value.evars;
						sntree[node].value.nodes.insert(sntree[node].value.nodes.end(),sntree[it].value.nodes.begin(),sntree[it].value.nodes.end());
						evars+=sntree[it].value.evars;
						int tmp=sntree[it].right_sibling;
						sntree.erase(sntree[it]);
						it=tmp;
						continue;
					}
					it=sntree[it].right_sibling;
				}
				if(sntree[node].right_sibling!=sntree.invalidPos) {
					node=sntree[node].right_sibling;
					down=true;
				}
			}
			else
				node=sntree[node].parent;
		}
		if(node==0||node==sntree.invalidPos)
			break;
	}
	return sntree;
}
}
}
}
#endif
