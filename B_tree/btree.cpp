#define _CRT_SECURE_NO_WARNINGS
#include "btree.h"
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

void BTProcessable::fromBytes(const std::vector<uint8_t> &obuf, size_t firstBytePos)
{
}
void BTProcessable::toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos)
{
}

BTProcessable &BTProcessable::operator=(const BTProcessable &other)
{
	return *this;
}

BTProcessable &BTProcessable::operator=(BTProcessable &&other)
{
	return *this;
}

BTProcessable *BTProcessable::createNew()
{
	return new BTProcessable();
}

bool BTProcessable::less(const BTProcessable &w2) const
{
	return false;
}

size_t BTProcessable::sizeInBytes()
{
	return 0;
}

bool operator<(const BTProcessable &w1, const BTProcessable &w2)
{
	return w1.less(w2);
}

bool operator==(const BTProcessable &w1, const BTProcessable &w2)
{
	return (!(w1 < w2 || w2 < w1));
}

bool operator<=(const BTProcessable &w1, const BTProcessable &w2)
{
	return (w1 < w2 || w2 == w1);
}

std::function<bool(const BTProcessable *, const BTProcessable *)> BTProcessablePtrComp = [](const BTProcessable *el1, const BTProcessable *el2)
{
	return (*el1 < *el2);
};

/////////////////////////////////////////////////////////////////

std::function<bool(const levelMarkersTreeNode &, const levelMarkersTreeNode &)> levelMarkersTreeNodeComp = [](const levelMarkersTreeNode &n1, const levelMarkersTreeNode &n2) -> bool
{
	return ((n1.levelMarker < n2.levelMarker) || (n1.levelMarker == n2.levelMarker && n1.id < n2.id));
};

/////////////////////////////////////////////////////////////////

NodeCache::NodeCache(B_Tree *ptrToTree, size_t maxNumOfPages) : _ptrToTree{ptrToTree}, _maxNumOfPages{maxNumOfPages}, _levelMarkersTree(levelMarkersTreeNodeComp) {}

bool NodeCache::processPretendentToCache(Node &node, bool getSave)
{
	if (_buffer.size() != _maxNumOfPages)
	{
		_buffer.insert({node._id, node});
		_levelMarkersTree.insert({node._levelMarker, node._id});

		return true;
	}

	auto minimum = *(_levelMarkersTree.begin());

	if (getSave == true)
	{
		if (minimum.levelMarker <= node._levelMarker)
		{
			_levelMarkersTree.erase(minimum);
			_levelMarkersTree.insert({node._levelMarker, node._id});

			writeNodeToDisk(_buffer.find(minimum.id)->second);

			_buffer.erase(minimum.id);
			_buffer.insert({node._id, node});

			return true;
		}
	}
	else
	{
		if (minimum.levelMarker <= node._levelMarker)
		{
			_levelMarkersTree.erase(minimum);
			_levelMarkersTree.insert({node._levelMarker, node._id});

			writeNodeToDisk(_buffer.find(minimum.id)->second);

			_buffer.erase(minimum.id);
			_buffer.insert({node._id, node});

			return true;
		}
	}

	return false;
}

void NodeCache::readNodeFromDisk(uint32_t id, Node &dst)
{
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	dst._id = id; // id

	std::vector<uint8_t> nodeBytes;
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	t1 = std::chrono::system_clock::now();
	_ptrToTree->nodesFile.seekg(id * _ptrToTree->_fileSizeInBytes);
	_ptrToTree->nodesFile.read((char *)&(nodeBytes[0]), nodeBytes.size());
	t2 = std::chrono::system_clock::now();

	////////////////////////////////

	size_t numOfPassedBytes = 0;

	dst.isLeaf = *((bool *)&(nodeBytes[numOfPassedBytes])); // isLeaf
	numOfPassedBytes += sizeof(bool);

	dst._numOfCurrentStoredObjects = *((uint16_t *)&(nodeBytes[numOfPassedBytes])); // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	dst._levelMarker = *((uint32_t *)&(nodeBytes[numOfPassedBytes])); // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	// for (size_t i = 0; i < dst.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	dst.childrenNodesIds[i] = *((uint32_t*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(uint32_t);
	// }
	std::memcpy((void *)&(dst.childrenNodesIds[0]), (void *)&(nodeBytes[numOfPassedBytes]), dst.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * dst.childrenNodesIds.size();

	for (size_t i = 0; i < dst.nodesValPtrs.size(); i++) // nodesVals
	{

		// dst.nodesValPtrs[i] = *((T *)&(nodeBytes[numOfPassedBytes]));
		dst.nodesValPtrs[i] = _ptrToTree->_initPtr->createNew();
		dst.nodesValPtrs[i]->fromBytes(nodeBytes, numOfPassedBytes);
		numOfPassedBytes += dst.nodesValPtrs[i]->sizeInBytes();
	}

	// std::memcpy((void *)&(dst.nodesVals[0]), (void *)&(nodeBytes[numOfPassedBytes]), dst.nodesVals.size() * sizeof(dst.nodesVals[0]));
	// numOfPassedBytes += sizeof(dst.nodesVals[0]) * dst.nodesVals.size();

	_ptrToTree->_rdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_rnum++;
}

void NodeCache::writeNodeToDisk(Node &src)
{
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	size_t numOfPassedBytes = 0;
	std::vector<uint8_t> nodeBytes;
	//	nodeBytes.assign((sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + sizeof(src.nodesVals[0]) * _ptrToTree->maxNumOfObjectsInNode + sizeof(src.childrenNodesIds[0]) * (_ptrToTree->maxNumOfObjectsInNode + 1)), 0);
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	*((bool *)&(nodeBytes[numOfPassedBytes])) = src.isLeaf; // isLeaf
	numOfPassedBytes += sizeof(bool);

	*((uint16_t *)&(nodeBytes[numOfPassedBytes])) = src._numOfCurrentStoredObjects; // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	*((uint32_t *)&(nodeBytes[numOfPassedBytes])) = src._levelMarker; // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	// for (size_t i = 0; i < src.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	*((uint32_t*)&(nodeBytes[numOfPassedBytes])) = src.childrenNodesIds[i];
	//	numOfPassedBytes += sizeof(uint32_t);
	// }
	std::memcpy((void *)&(nodeBytes[numOfPassedBytes]), (void *)&(src.childrenNodesIds[0]), src.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * src.childrenNodesIds.size();

	for (size_t i = 0; i < src.nodesValPtrs.size(); i++) // nodesVals
	{
		src.nodesValPtrs[i]->toBytes(nodeBytes, numOfPassedBytes);
		numOfPassedBytes += src.nodesValPtrs[i]->sizeInBytes();
	}
	// std::memcpy((void *)&(nodeBytes[numOfPassedBytes]), (void *)&(src.nodesVals[0]), src.nodesVals.size() * sizeof(src.nodesVals[0]));
	// numOfPassedBytes += sizeof(T) * src.nodesVals.size();

	////////////////////////////////

	t1 = std::chrono::system_clock::now();
	_ptrToTree->nodesFile.seekp(src._id * _ptrToTree->_fileSizeInBytes);
	_ptrToTree->nodesFile.write((char *)&(nodeBytes[0]), nodeBytes.size());
	t2 = std::chrono::system_clock::now();

	_ptrToTree->_wdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_wnum++;
}

void NodeCache::getNode(uint32_t id, Node &dst)
{
	if (_maxNumOfPages == 0)
	{
		readNodeFromDisk(id, dst);
		return;
	}

	auto iter = _buffer.find(id);

	if (iter != _buffer.end())
	{
		dst = iter->second;
	}
	else
	{
		readNodeFromDisk(id, dst);
		processPretendentToCache(dst, true);
	}
}

void NodeCache::saveNode(Node &src)
{
	if (_maxNumOfPages == 0)
	{
		writeNodeToDisk(src);
		return;
	}

	auto iter = _buffer.find(src._id);

	if (iter != _buffer.end())
	{
		iter->second = src; // levelMarker must not be changed
	}
	else
	{
		if (processPretendentToCache(src, false) == false)
			writeNodeToDisk(src);
	}
}

void NodeCache::freeNode(Node &node)
{
	// remove from _freqTree and from buffer
	if (_buffer.find(node._id) != _buffer.end())
	{
		_levelMarkersTree.erase({node._levelMarker, node._id});
		_buffer.erase(node._id);
	}
}

/////////////////////////////////////////////////////////////////

B_Tree::B_Tree(BTProcessable *initPtr, size_t objectSizeInBytes, bool recover, bool keepOnDiskAfterDestruction, std::string _pathToNodes, size_t diskPageSizeInBytes, size_t maxNumOfPagesInCache) : maxNumOfObjectsInNode{uint16_t((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (objectSizeInBytes + sizeof(uint32_t)) - ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (objectSizeInBytes + sizeof(uint32_t)) + 1) % 2)},
																																																	 _root(initPtr, (diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (objectSizeInBytes + sizeof(uint32_t)) - ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (objectSizeInBytes + sizeof(uint32_t)) + 1) % 2, 0),
																																																	 pathToNodes{_pathToNodes},
																																																	 _initPtr{initPtr},
																																																	 _keepOnDiskAfterDestruction{keepOnDiskAfterDestruction},
																																																	 _cache(this, maxNumOfPagesInCache),
																																																	 _rdur(0),
																																																	 _wdur(0),
																																																	 _rnum(0),
																																																	 _wnum(0)
{
	t = (maxNumOfObjectsInNode + 1) / 2;
	minNumOfObjectsInNode = t - 1;
	largestExistingId = 0;
	_fileSizeInBytes = sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + objectSizeInBytes * maxNumOfObjectsInNode + sizeof(uint32_t) * (maxNumOfObjectsInNode + 1);

	//////
	if (pathToNodes.empty())
	{
		pathToNodes = "./";
	}

	nodesFileName = pathToNodes + std::string("file.btr");
	std::ofstream(nodesFileName.c_str()).close();
	nodesFile.open(nodesFileName.c_str(), std::ios::binary | std::ios::out | std::ios::in);
	if (!nodesFile)
	{
		std::cerr << "\"file.btr\" could not be opened!" << std::endl;
		exit(1);
	}
	//////

	// std::cout << "\n<< " << maxNumOfObjectsInNode << " objects in node>>\n";
	createTree();
}

void B_Tree::getNode(uint32_t id, Node &dst)
{
	if (id == _root._id)
	{
		dst = _root;
	}
	else
	{
		_cache.getNode(id, dst);
	}
}

void B_Tree::saveNode(Node &src)
{
	if (src._id == _root._id)
		_root = src;
	else
	{
		_cache.saveNode(src);
	}
}

void B_Tree::freeNode(Node &node)
{
	if (node._id != _root._id)
	{
		_cache.freeNode(node);
	}
}

void B_Tree::makeRootUpper()
{
	Node newRoot(_initPtr, maxNumOfObjectsInNode, _root._levelMarker + 1);
	initNode(newRoot);

	newRoot.isLeaf = false;
	newRoot._numOfCurrentStoredObjects = 0;
	newRoot.childrenNodesIds[0] = _root._id;

	Node prevRoot = _root;

	_root = newRoot;

	saveNode(prevRoot);
}

void B_Tree::makeRootLower(Node &newRoot)
{
	freeNode(newRoot);
	_root = newRoot;
}

void B_Tree::createTree()
{
	initNode(_root);
	_root.isLeaf = true;
	_root._numOfCurrentStoredObjects = 0;
}

void B_Tree::search(Node &root, BTProcessable &val, Node &dst_node, uint16_t &indexWhereValueFound)
{
	uint16_t i = 1;

	// while ((i <= root._numOfCurrentStoredObjects) && (val > root.nodesVals[i - 1]))
	//	i++;
	BTProcessable *valCopyPtr = val.createNew();
	*valCopyPtr = val;
	auto it = std::upper_bound(root.nodesValPtrs.begin(), root.nodesValPtrs.begin() + root._numOfCurrentStoredObjects, valCopyPtr, BTProcessablePtrComp);
	if (it != root.nodesValPtrs.begin() && *(*(it - 1)) == val)
		it = it - 1;
	i += (it - root.nodesValPtrs.begin());

	if ((i <= root._numOfCurrentStoredObjects) && (val == *root.nodesValPtrs[i - 1]))
	{
		dst_node = root;
		indexWhereValueFound = i;
	}
	else if (root.isLeaf)
		indexWhereValueFound = 0;
	else
	{
		Node newRoot(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1);
		getNode(root.childrenNodesIds[i - 1], newRoot);
		search(newRoot, val, dst_node, indexWhereValueFound);
	}
}

void B_Tree::splitChild(Node &node, uint16_t index)
{
	Node newRight(_initPtr, maxNumOfObjectsInNode, node._levelMarker - 1);
	initNode(newRight);

	Node splitingChild(_initPtr, maxNumOfObjectsInNode, node._levelMarker - 1);
	getNode(node.childrenNodesIds[index - 1], splitingChild);

	newRight.isLeaf = splitingChild.isLeaf;
	newRight._numOfCurrentStoredObjects = minNumOfObjectsInNode;

	for (uint16_t j = 1; j <= minNumOfObjectsInNode; j++)
		newRight.nodesValPtrs[j - 1] = splitingChild.nodesValPtrs[t + j - 1]; // min num of objs in left, one up, min num of objs in right

	if (splitingChild.isLeaf == false)
	{
		for (uint16_t j = 1; j <= minNumOfObjectsInNode + 1; j++)
			newRight.childrenNodesIds[j - 1] = splitingChild.childrenNodesIds[t + j - 1]; // min num of childs in left, min num of childs in right
	}

	splitingChild._numOfCurrentStoredObjects = minNumOfObjectsInNode;

	for (uint16_t j = node._numOfCurrentStoredObjects + 1; j >= index + 1; j--)
		node.childrenNodesIds[j + 1 - 1] = node.childrenNodesIds[j - 1];
	node.childrenNodesIds[index + 1 - 1] = newRight._id;

	for (uint16_t j = node._numOfCurrentStoredObjects; j >= index; j--)
		node.nodesValPtrs[j + 1 - 1] = node.nodesValPtrs[j - 1];
	node.nodesValPtrs[index - 1] = splitingChild.nodesValPtrs[t - 1];

	node._numOfCurrentStoredObjects++;

	saveNode(splitingChild);
	saveNode(newRight);
	saveNode(node);
}

void B_Tree::insertIntoNonFull(Node &node, BTProcessable &val)
{
	auto i = node._numOfCurrentStoredObjects;

	if (node.isLeaf)
	{
		while (i >= 1 && val < *(node.nodesValPtrs[i - 1]))
		{
			node.nodesValPtrs[i + 1 - 1] = node.nodesValPtrs[i - 1];
			i--;
		}

		BTProcessable *valCopyPtr{val.createNew()};
		*valCopyPtr = val;
		node.nodesValPtrs[i + 1 - 1] = valCopyPtr;
		node._numOfCurrentStoredObjects++;

		saveNode(node);
	}
	else
	{
		while (i >= 1 && val < (*node.nodesValPtrs[i - 1]))
			i--;
		i++;

		Node nextNode(_initPtr, maxNumOfObjectsInNode, node._levelMarker - 1);
		getNode(node.childrenNodesIds[i - 1], nextNode);

		if (nextNode._numOfCurrentStoredObjects == maxNumOfObjectsInNode)
		{
			splitChild(node, i);

			if (*node.nodesValPtrs[i - 1] < val)
				i++;
		}

		Node nextNode1(_initPtr, maxNumOfObjectsInNode, node._levelMarker - 1);
		getNode(node.childrenNodesIds[i - 1], nextNode1);

		insertIntoNonFull(nextNode1, val);
	}
}

void B_Tree::insertKey(BTProcessable &val)
{
	if (_root._numOfCurrentStoredObjects == maxNumOfObjectsInNode)
	{
		makeRootUpper();
		splitChild(_root, 1); // previos root will be saved here
	}

	insertIntoNonFull(_root, val);
}

void B_Tree::deleteKey(Node &root, BTProcessable &val)
{
	uint16_t foundIndex = 0;

	for (uint16_t i = 0; i < root._numOfCurrentStoredObjects; i++)
	{
		if ((*root.nodesValPtrs[root._numOfCurrentStoredObjects - 1 - i]) <= val)
		{
			foundIndex = root._numOfCurrentStoredObjects - 1 - i + 1;
			break;
		}
	}

	if ((foundIndex != 0) && (*root.nodesValPtrs[foundIndex - 1] == val)) // 1, 2
	{
		foundIndex--;

		if (root.isLeaf == true) // 1
		{
			for (uint16_t i = foundIndex; i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				root.nodesValPtrs[i] = root.nodesValPtrs[i + 1];

			root._numOfCurrentStoredObjects--;

			saveNode(root);
		}
		else // 2
		{
			Node leftNeighbourNode(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1), rightNeighbourNode(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1);

			getNode(root.childrenNodesIds[foundIndex], leftNeighbourNode);		// left
			getNode(root.childrenNodesIds[foundIndex + 1], rightNeighbourNode); // right

			if (leftNeighbourNode._numOfCurrentStoredObjects >= t) // a
			{
				while (leftNeighbourNode.isLeaf != true)
					getNode(leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects], leftNeighbourNode);

				BTProcessable *valToMove = (*leftNeighbourNode.nodesValPtrs[leftNeighbourNode._numOfCurrentStoredObjects - 1]).createNew();
				*valToMove = *leftNeighbourNode.nodesValPtrs[leftNeighbourNode._numOfCurrentStoredObjects - 1];
				// BTProcessable valToMove2 = *leftNeighbourNode.nodesValPtrs[leftNeighbourNode._numOfCurrentStoredObjects - 1];
				root.nodesValPtrs[foundIndex] = valToMove; // copying - we need instance further

				saveNode(root);

				getNode(root.childrenNodesIds[foundIndex], leftNeighbourNode);
				deleteKey(leftNeighbourNode, *valToMove);
			}
			else if (rightNeighbourNode._numOfCurrentStoredObjects >= t) // b
			{
				while (rightNeighbourNode.isLeaf != true)
					getNode(rightNeighbourNode.childrenNodesIds[0], rightNeighbourNode);

				BTProcessable *valToMove = (*rightNeighbourNode.nodesValPtrs[0]).createNew();
				*valToMove = *rightNeighbourNode.nodesValPtrs[0];
				// BTProcessable valToMove2 = *rightNeighbourNode.nodesValPtrs[0];

				root.nodesValPtrs[foundIndex] = valToMove; // copying - we need instance further

				saveNode(root);

				getNode(root.childrenNodesIds[foundIndex + 1], rightNeighbourNode);
				deleteKey(rightNeighbourNode, *valToMove);
			}
			else // c
			{
				BTProcessable *valCopy = val.createNew();
				*valCopy = val;
				leftNeighbourNode.nodesValPtrs[leftNeighbourNode._numOfCurrentStoredObjects] = valCopy; // copying - we need instance further
				leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects + 1] = rightNeighbourNode.childrenNodesIds[0];

				for (size_t i = 0; i < rightNeighbourNode._numOfCurrentStoredObjects; i++)
				{
					leftNeighbourNode.nodesValPtrs[leftNeighbourNode._numOfCurrentStoredObjects + 1 + i] = rightNeighbourNode.nodesValPtrs[i];
					leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects + 1 + 1 + i] = rightNeighbourNode.childrenNodesIds[i + 1];
				}

				leftNeighbourNode._numOfCurrentStoredObjects += (1 + rightNeighbourNode._numOfCurrentStoredObjects);

				for (uint16_t i = foundIndex; i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				{
					root.nodesValPtrs[i] = root.nodesValPtrs[i + 1];
					root.childrenNodesIds[i + 1] = root.childrenNodesIds[i + 1 + 1];
				}

				root._numOfCurrentStoredObjects--;

				saveNode(root);

				if (root._id == _root._id)
				{
					if (root._numOfCurrentStoredObjects == 0)
						makeRootLower(leftNeighbourNode);
				}

				saveNode(leftNeighbourNode);

				freeNode(rightNeighbourNode);

				deleteKey(leftNeighbourNode, val);
			}
		}
	}
	else // 3
	{
		Node mayBeContain(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1);

		uint16_t mayBeContainRefNumber = (foundIndex != 0) ? (foundIndex - 1 + 1) : 0; // = foundIndex

		getNode(root.childrenNodesIds[mayBeContainRefNumber], mayBeContain);

		// foundIndex--;

		if (mayBeContain._numOfCurrentStoredObjects == minNumOfObjectsInNode)
		{
			Node neigbourOfMayBeContain(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1);

			bool foundNeigbourOfMayBeContainWithTKeys = false;
			uint16_t neigbourOfMayBeContainRefNumber = 0;

			if (mayBeContainRefNumber > 0)
			{
				neigbourOfMayBeContainRefNumber = mayBeContainRefNumber - 1; // left
				getNode(root.childrenNodesIds[neigbourOfMayBeContainRefNumber], neigbourOfMayBeContain);

				if (neigbourOfMayBeContain._numOfCurrentStoredObjects >= t)
					foundNeigbourOfMayBeContainWithTKeys = true; // a
			}

			if ((mayBeContainRefNumber < root._numOfCurrentStoredObjects) && (foundNeigbourOfMayBeContainWithTKeys == false))
			{
				neigbourOfMayBeContainRefNumber = mayBeContainRefNumber + 1; // right
				getNode(root.childrenNodesIds[neigbourOfMayBeContainRefNumber], neigbourOfMayBeContain);

				if (neigbourOfMayBeContain._numOfCurrentStoredObjects >= t)
					foundNeigbourOfMayBeContainWithTKeys = true; // a
			}

			if (foundNeigbourOfMayBeContainWithTKeys == true) // a
			{
				if (neigbourOfMayBeContainRefNumber > mayBeContainRefNumber)
				{
					mayBeContain.nodesValPtrs[mayBeContain._numOfCurrentStoredObjects] = root.nodesValPtrs[mayBeContainRefNumber];
					mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1] = neigbourOfMayBeContain.childrenNodesIds[0];
					mayBeContain._numOfCurrentStoredObjects++;

					root.nodesValPtrs[mayBeContainRefNumber] = neigbourOfMayBeContain.nodesValPtrs[0];

					neigbourOfMayBeContain._numOfCurrentStoredObjects--;

					for (size_t i = 0; i < neigbourOfMayBeContain._numOfCurrentStoredObjects; i++)
					{
						neigbourOfMayBeContain.nodesValPtrs[i] = neigbourOfMayBeContain.nodesValPtrs[i + 1];
						neigbourOfMayBeContain.childrenNodesIds[i] = neigbourOfMayBeContain.childrenNodesIds[i + 1];
					}

					neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects] = neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects + 1];
				}
				else
				{
					for (size_t i = 0; i < mayBeContain._numOfCurrentStoredObjects; i++)
					{
						mayBeContain.nodesValPtrs[mayBeContain._numOfCurrentStoredObjects - 1 - i + 1] = mayBeContain.nodesValPtrs[mayBeContain._numOfCurrentStoredObjects - 1 - i];
						mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1 - 1 - i + 1] = mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1 - 1 - i];
					}
					mayBeContain.childrenNodesIds[1] = mayBeContain.childrenNodesIds[0];

					mayBeContain.nodesValPtrs[0] = root.nodesValPtrs[neigbourOfMayBeContainRefNumber];
					mayBeContain.childrenNodesIds[0] = neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects];
					mayBeContain._numOfCurrentStoredObjects++;

					root.nodesValPtrs[neigbourOfMayBeContainRefNumber] = neigbourOfMayBeContain.nodesValPtrs[neigbourOfMayBeContain._numOfCurrentStoredObjects - 1];

					neigbourOfMayBeContain._numOfCurrentStoredObjects--;
				}

				saveNode(root);
				saveNode(neigbourOfMayBeContain);
				saveNode(mayBeContain);

				deleteKey(mayBeContain, val);
			}
			else // b
			{
				Node *left;
				Node *right;

				if (neigbourOfMayBeContainRefNumber > mayBeContainRefNumber)
				{
					left = &mayBeContain;
					right = &neigbourOfMayBeContain;
				}
				else
				{
					left = &neigbourOfMayBeContain;
					right = &mayBeContain;
				}

				left->nodesValPtrs[left->_numOfCurrentStoredObjects] = root.nodesValPtrs[std::min(mayBeContainRefNumber, neigbourOfMayBeContainRefNumber)];
				left->childrenNodesIds[left->_numOfCurrentStoredObjects + 1] = right->childrenNodesIds[0];

				for (size_t i = 0; i < right->_numOfCurrentStoredObjects; i++)
				{
					left->nodesValPtrs[left->_numOfCurrentStoredObjects + 1 + i] = right->nodesValPtrs[i];
					left->childrenNodesIds[left->_numOfCurrentStoredObjects + 1 + 1 + i] = right->childrenNodesIds[i + 1];
				}

				left->_numOfCurrentStoredObjects += (1 + right->_numOfCurrentStoredObjects);

				for (uint16_t i = std::min(mayBeContainRefNumber, neigbourOfMayBeContainRefNumber); i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				{
					root.nodesValPtrs[i] = root.nodesValPtrs[i + 1];
					root.childrenNodesIds[i + 1] = root.childrenNodesIds[i + 1 + 1];
				}

				root._numOfCurrentStoredObjects--;

				saveNode(root);

				if (root._id == _root._id)
				{
					if (root._numOfCurrentStoredObjects == 0)
					{
						if (neigbourOfMayBeContainRefNumber == 0)
						{
							makeRootLower(neigbourOfMayBeContain);
						}
						else if (mayBeContainRefNumber == 0)
						{
							makeRootLower(mayBeContain);
						}
						else
						{
							std::cout << "\n========\nIMPOSSIBLE!!!\n========\n";
							Node newRoot(_initPtr, maxNumOfObjectsInNode, root._levelMarker - 1);
							getNode(_root.childrenNodesIds[0], newRoot);

							makeRootLower(newRoot);
						}
					}
				}

				if (neigbourOfMayBeContainRefNumber > mayBeContainRefNumber)
				{
					saveNode(mayBeContain);
					freeNode(neigbourOfMayBeContain);

					deleteKey(mayBeContain, val);
				}
				else
				{
					saveNode(neigbourOfMayBeContain);
					freeNode(mayBeContain);

					deleteKey(neigbourOfMayBeContain, val);
				}
			}
		}
		else // mayBeContain has t objects
		{
			deleteKey(mayBeContain, val);
		}
	}
}

B_Tree::~B_Tree()
{
	nodesFile.close();

	if (_keepOnDiskAfterDestruction == false)
		std::filesystem::remove(nodesFileName);
}