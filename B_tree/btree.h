#ifndef BTREE
#define BTREE

// #define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <thread>
#include <set>
#include <cstring>
#include <algorithm>
#include <stdio.h>
#include <functional>

/////////////////////////////////////////////////////////////////

template <class T>
class Node
{
public:
	Node(uint16_t numOfObjectsInNode, uint32_t levelMarker) : _levelMarker{levelMarker}, _id{0}, isLeaf{false}, _numOfCurrentStoredObjects{0}
	{
		childrenNodesIds.assign(numOfObjectsInNode + 1, 0);
		nodesVals.assign(numOfObjectsInNode, T());
	}

	uint32_t _id;
	bool isLeaf;
	uint16_t _numOfCurrentStoredObjects;
	uint32_t _levelMarker;

	std::vector<uint32_t> childrenNodesIds;
	std::vector<T> nodesVals;

	// in memory:									  - id - | - is leaf - | - num OfCurrent Stored Objects - | - level marker - | - children Node's Ids - | - nodes values -
	// in disk: filename = <id>.btrnd; filecontent =	       - is leaf - | - num OfCurrent Stored Objects - | - level marker - | - children Node's Ids - | - nodes values -
};

/////////////////////////////////////////////////////////////////

struct levelMarkersTreeNode
{
	uint32_t levelMarker;
	uint32_t id;
};

/////////////////////////////////////////////////////////////////

template <class T>
class B_Tree;

template <class T>
class NodeCache
{
public:
	size_t _maxNumOfPages;
	std::unordered_map<uint32_t, Node<T>> _buffer;
	std::set<levelMarkersTreeNode, std::function<bool(const levelMarkersTreeNode &, const levelMarkersTreeNode &)>> _levelMarkersTree;
	B_Tree<T> *_ptrToTree;

	NodeCache(B_Tree<T> *ptrToTree, size_t maxNumOfPages);

	bool processPretendentToCache(Node<T> &node, bool getSave);
	void getNode(uint32_t id, Node<T> &dst);
	void saveNode(Node<T> &src);
	void freeNode(Node<T> &node);

	void writeNodeToDisk(Node<T> &src);
	void readNodeFromDisk(uint32_t id, Node<T> &dst);
};

/////////////////////////////////////////////////////////////////

template <class T>
class B_Tree
{
public:
	std::chrono::microseconds _rdur;
	std::chrono::microseconds _wdur;
	uint64_t _rnum;
	uint64_t _wnum;

	uint32_t largestExistingId;
	const uint16_t maxNumOfObjectsInNode; // = 2 * t - 1
	uint16_t minNumOfObjectsInNode;		  // = t - 1
	uint16_t t;
	std::string pathToNodes;
	std::string nodesFileName;
	bool _keepOnDiskAfterDestruction;
	int64_t _fileSizeInBytes;
	std::fstream nodesFile;

	NodeCache<T> _cache;
	Node<T> _root;

	B_Tree(bool recover = false, bool keepOnDiskAfterDestruction = false, std::string _pathToNodes = "", size_t diskPageSizeInBytes = 4096, size_t maxNumOfPagesInCache = 4);
	~B_Tree();

	void getNode(uint32_t id, Node<T> &dst);
	void saveNode(Node<T> &src);
	void initNode(Node<T> &node)
	{
		node._id = largestExistingId + 1;
		largestExistingId += 1;
	}
	void freeNode(Node<T> &node);
	void makeRootUpper();
	void makeRootLower(Node<T> &newRoot);
	void createTree();
	void search(Node<T> &root, T &val, Node<T> &dst_node, uint16_t &indexWhereValueFound);
	void splitChild(Node<T> &node, uint16_t index);
	void insertKey(T &val);
	void insertIntoNonFull(Node<T> &node, T &val);
	void deleteKey(Node<T> &root, T &val);
};

#endif