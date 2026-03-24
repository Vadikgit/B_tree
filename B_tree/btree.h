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
class BTProcessable
{
public:
	virtual void fromBytes(const std::vector<uint8_t> &ibuf, size_t firstBytePos);
	virtual void toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos);
	virtual bool less(const BTProcessable &w2) const;
	virtual size_t sizeInBytes();

	BTProcessable() = default;
	BTProcessable(BTProcessable &&other) = default;
	BTProcessable(const BTProcessable &other) = default;
	virtual ~BTProcessable() {};

	virtual BTProcessable *createNew();

	virtual BTProcessable &operator=(const BTProcessable &other);
	virtual BTProcessable &operator=(BTProcessable &&other);
};

bool operator<(const BTProcessable &w1, const BTProcessable &w2);
bool operator==(const BTProcessable &w1, const BTProcessable &w2);
bool operator<=(const BTProcessable &w1, const BTProcessable &w2);

class Node
{
public:
	Node(BTProcessable *initPtr, uint16_t numOfObjectsInNode, uint32_t levelMarker) : _levelMarker{levelMarker}, _id{0}, isLeaf{false}, _numOfCurrentStoredObjects{0}
	{
		childrenNodesIds.assign(numOfObjectsInNode + 1, 0);

		for (size_t i = 0; i < numOfObjectsInNode; i++)
		{
			nodesValPtrs.push_back(std::unique_ptr<BTProcessable>{initPtr->createNew()});
		}
	}

	Node(const Node &other) : _id{other._id}, isLeaf{other.isLeaf}, _numOfCurrentStoredObjects{other._numOfCurrentStoredObjects}, _levelMarker{other._levelMarker}, childrenNodesIds{other.childrenNodesIds}
	{
		for (size_t i = 0; i < other.nodesValPtrs.size(); i++)
		{
			nodesValPtrs.push_back(std::unique_ptr<BTProcessable>{other.nodesValPtrs[i]->createNew()});
			*nodesValPtrs[i] = *other.nodesValPtrs[i];
		}
	}

	Node &operator=(const Node &other)
	{
		if (&other != this)
		{
			_id = other._id;
			isLeaf = other.isLeaf;
			_numOfCurrentStoredObjects = other._numOfCurrentStoredObjects;
			_levelMarker = other._levelMarker;
			childrenNodesIds = other.childrenNodesIds;

			for (size_t i = 0; i < other.nodesValPtrs.size(); i++)
			{
				nodesValPtrs[i].reset(other.nodesValPtrs[i]->createNew());
				*nodesValPtrs[i] = *other.nodesValPtrs[i];
			}
		}

		return *this;
	}

	Node &operator=(Node &&other)
	{
		_id = other._id;
		isLeaf = other.isLeaf;
		_numOfCurrentStoredObjects = other._numOfCurrentStoredObjects;
		_levelMarker = other._levelMarker;
		childrenNodesIds = std::move(other.childrenNodesIds);
		nodesValPtrs = std::move(other.nodesValPtrs);

		return *this;
	}

	Node(Node &&other) : _id{other._id}, isLeaf{other.isLeaf}, _numOfCurrentStoredObjects{other._numOfCurrentStoredObjects}, _levelMarker{other._levelMarker}, childrenNodesIds{std::move(other.childrenNodesIds)}, nodesValPtrs{std::move(other.nodesValPtrs)}
	{
	}

	uint32_t _id;
	bool isLeaf;
	uint16_t _numOfCurrentStoredObjects;
	uint32_t _levelMarker;

	std::vector<uint32_t> childrenNodesIds;
	std::vector<std::unique_ptr<BTProcessable>> nodesValPtrs;

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

class B_Tree;

class NodeCache
{
public:
	size_t _maxNumOfPages;
	std::unordered_map<uint32_t, Node> _buffer;
	std::set<levelMarkersTreeNode, std::function<bool(const levelMarkersTreeNode &, const levelMarkersTreeNode &)>> _levelMarkersTree;
	B_Tree *_ptrToTree;

	NodeCache(B_Tree *ptrToTree, size_t maxNumOfPages);

	bool processPretendentToCache(Node &node, bool getSave);
	void getNode(uint32_t id, Node &dst);
	void saveNode(Node &src);
	void freeNode(Node &node);

	void writeNodeToDisk(Node &src);
	void readNodeFromDisk(uint32_t id, Node &dst);
};

/////////////////////////////////////////////////////////////////

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

	NodeCache _cache;
	Node _root;
	BTProcessable *_initPtr;

	B_Tree(BTProcessable *initPtr, size_t objectSizeInBytes, bool recover = false, bool keepOnDiskAfterDestruction = false, std::string _pathToNodes = "", size_t diskPageSizeInBytes = 4096, size_t maxNumOfPagesInCache = 4);
	~B_Tree();

	void getNode(uint32_t id, Node &dst);
	void saveNode(Node &src);
	void initNode(Node &node)
	{
		node._id = largestExistingId + 1;
		largestExistingId += 1;
	}
	void freeNode(Node &node);
	void makeRootUpper();
	void makeRootLower(Node &newRoot);
	void createTree();
	void search(Node &root, BTProcessable &val, Node &dst_node, uint16_t &indexWhereValueFound);
	void splitChild(Node &node, uint16_t index);
	void insertKey(BTProcessable &val);
	void insertIntoNonFull(Node &node, BTProcessable &val);
	void deleteKey(Node &root, BTProcessable &val);
};

#endif