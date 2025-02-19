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

/////////////////////////////////////////////////////////////////

template class B_Tree<uint64_t>;

/////////////////////////////////////////////////////////////////

std::function<bool(const levelMarkersTreeNode &, const levelMarkersTreeNode &)> levelMarkersTreeNodeComp = [](const levelMarkersTreeNode &n1, const levelMarkersTreeNode &n2) -> bool
{
	return ((n1.levelMarker < n2.levelMarker) || (n1.levelMarker == n2.levelMarker && n1.id < n2.id));
};

/////////////////////////////////////////////////////////////////

template <class T>
NodeCache<T>::NodeCache(B_Tree<T> *ptrToTree, size_t maxNumOfPages) : _ptrToTree{ptrToTree}, _maxNumOfPages{maxNumOfPages}, _levelMarkersTree(levelMarkersTreeNodeComp) {}

template <class T>
bool NodeCache<T>::processPretendentToCache(Node<T> &node, bool getSave)
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

template <class T>
void NodeCache<T>::readNodeFromDisk(uint32_t id, Node<T> &dst)
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

	// for (size_t i = 0; i < dst.nodesVals.size(); i++) // nodesVals
	//{
	//	dst.nodesVals[i] = *((T*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(T);
	// }
	std::memcpy((void *)&(dst.nodesVals[0]), (void *)&(nodeBytes[numOfPassedBytes]), dst.nodesVals.size() * sizeof(dst.nodesVals[0]));
	numOfPassedBytes += sizeof(dst.nodesVals[0]) * dst.nodesVals.size();

	_ptrToTree->_rdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_rnum++;
}

template <class T>
void NodeCache<T>::writeNodeToDisk(Node<T> &src)
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

	// for (size_t i = 0; i < src.nodesVals.size(); i++) // nodesVals
	//{
	//	*((T*)&(nodeBytes[numOfPassedBytes])) = src.nodesVals[i];
	//	numOfPassedBytes += sizeof(T);
	// }
	std::memcpy((void *)&(nodeBytes[numOfPassedBytes]), (void *)&(src.nodesVals[0]), src.nodesVals.size() * sizeof(src.nodesVals[0]));
	numOfPassedBytes += sizeof(T) * src.nodesVals.size();

	////////////////////////////////

	t1 = std::chrono::system_clock::now();
	_ptrToTree->nodesFile.seekp(src._id * _ptrToTree->_fileSizeInBytes);
	_ptrToTree->nodesFile.write((char *)&(nodeBytes[0]), nodeBytes.size());
	t2 = std::chrono::system_clock::now();

	_ptrToTree->_wdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_wnum++;
}

template <class T>
void NodeCache<T>::getNode(uint32_t id, Node<T> &dst)
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

template <class T>
void NodeCache<T>::saveNode(Node<T> &src)
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

template <class T>
void NodeCache<T>::freeNode(Node<T> &node)
{
	// remove from _freqTree and from buffer
	if (_buffer.find(node._id) != _buffer.end())
	{
		_levelMarkersTree.erase({node._levelMarker, node._id});
		_buffer.erase(node._id);
	}
}

/////////////////////////////////////////////////////////////////

template <class T>
B_Tree<T>::B_Tree(bool recover, bool keepOnDiskAfterDestruction, std::string _pathToNodes, size_t diskPageSizeInBytes, size_t maxNumOfPagesInCache) : maxNumOfObjectsInNode{uint16_t((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) - ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) + 1) % 2)},
																																					  _root((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) - ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) + 1) % 2, 0),
																																					  pathToNodes{_pathToNodes},
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
	_fileSizeInBytes = sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + sizeof(T) * maxNumOfObjectsInNode + sizeof(uint32_t) * (maxNumOfObjectsInNode + 1);

	//////
	nodesFileName = pathToNodes + std::string("/file.btr");
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

template <class T>
void B_Tree<T>::getNode(uint32_t id, Node<T> &dst)
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

template <class T>
void B_Tree<T>::saveNode(Node<T> &src)
{
	if (src._id == _root._id)
		_root = src;
	else
	{
		_cache.saveNode(src);
	}
}

template <class T>
void B_Tree<T>::freeNode(Node<T> &node)
{
	if (node._id != _root._id)
	{
		_cache.freeNode(node);
	}
}

template <class T>
void B_Tree<T>::makeRootUpper()
{
	Node<T> newRoot(maxNumOfObjectsInNode, _root._levelMarker + 1);
	initNode(newRoot);

	newRoot.isLeaf = false;
	newRoot._numOfCurrentStoredObjects = 0;
	newRoot.childrenNodesIds[0] = _root._id;

	Node<T> prevRoot = _root;

	_root = newRoot;

	saveNode(prevRoot);
}

template <class T>
void B_Tree<T>::makeRootLower(Node<T> &newRoot)
{
	freeNode(newRoot);
	_root = newRoot;
}

template <class T>
void B_Tree<T>::createTree()
{
	initNode(_root);
	_root.isLeaf = true;
	_root._numOfCurrentStoredObjects = 0;
}

template <class T>
void B_Tree<T>::search(Node<T> &root, T &val, Node<T> &dst_node, uint16_t &indexWhereValueFound)
{
	uint16_t i = 1;

	// while ((i <= root._numOfCurrentStoredObjects) && (val > root.nodesVals[i - 1]))
	//	i++;
	auto it = std::upper_bound(root.nodesVals.begin(), root.nodesVals.begin() + root._numOfCurrentStoredObjects, val);
	if (it != root.nodesVals.begin() && *(it - 1) == val)
		it = it - 1;
	i += (it - root.nodesVals.begin());

	if ((i <= root._numOfCurrentStoredObjects) && (val == root.nodesVals[i - 1]))
	{
		dst_node = root;
		indexWhereValueFound = i;
	}
	else if (root.isLeaf)
		indexWhereValueFound = 0;
	else
	{
		Node<T> newRoot(maxNumOfObjectsInNode, root._levelMarker - 1);
		getNode(root.childrenNodesIds[i - 1], newRoot);
		search(newRoot, val, dst_node, indexWhereValueFound);
	}
}

template <class T>
void B_Tree<T>::splitChild(Node<T> &node, uint16_t index)
{
	Node<T> newRight(maxNumOfObjectsInNode, node._levelMarker - 1);
	initNode(newRight);

	Node<T> splitingChild(maxNumOfObjectsInNode, node._levelMarker - 1);
	getNode(node.childrenNodesIds[index - 1], splitingChild);

	newRight.isLeaf = splitingChild.isLeaf;
	newRight._numOfCurrentStoredObjects = minNumOfObjectsInNode;

	for (uint16_t j = 1; j <= minNumOfObjectsInNode; j++)
		newRight.nodesVals[j - 1] = splitingChild.nodesVals[t + j - 1]; // min num of objs in left, one up, min num of objs in right

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
		node.nodesVals[j + 1 - 1] = node.nodesVals[j - 1];
	node.nodesVals[index - 1] = splitingChild.nodesVals[t - 1];

	node._numOfCurrentStoredObjects++;

	saveNode(splitingChild);
	saveNode(newRight);
	saveNode(node);
}

template <class T>
void B_Tree<T>::insertIntoNonFull(Node<T> &node, T &val)
{
	auto i = node._numOfCurrentStoredObjects;

	if (node.isLeaf)
	{
		while (i >= 1 && val < node.nodesVals[i - 1])
		{
			node.nodesVals[i - 1 + 1] = node.nodesVals[i - 1];
			i--;
		}

		node.nodesVals[i - 1 + 1] = val;
		node._numOfCurrentStoredObjects++;

		saveNode(node);
	}
	else
	{
		while (i >= 1 && val < node.nodesVals[i - 1])
			i--;
		i++;

		Node<T> nextNode(maxNumOfObjectsInNode, node._levelMarker - 1);
		getNode(node.childrenNodesIds[i - 1], nextNode);

		if (nextNode._numOfCurrentStoredObjects == maxNumOfObjectsInNode)
		{
			splitChild(node, i);

			if (val > node.nodesVals[i - 1])
				i++;
		}

		Node<T> nextNode1(maxNumOfObjectsInNode, node._levelMarker - 1);
		getNode(node.childrenNodesIds[i - 1], nextNode1);

		insertIntoNonFull(nextNode1, val);
	}
}

template <class T>
void B_Tree<T>::insertKey(T &val)
{
	if (_root._numOfCurrentStoredObjects == maxNumOfObjectsInNode)
	{
		makeRootUpper();
		splitChild(_root, 1); // previos root will be saved here
	}

	insertIntoNonFull(_root, val);
}

template <class T>
void B_Tree<T>::deleteKey(Node<T> &root, T &val)
{
	uint16_t foundIndex = 0;

	for (uint16_t i = 0; i < root._numOfCurrentStoredObjects; i++)
	{
		if (root.nodesVals[root._numOfCurrentStoredObjects - 1 - i] <= val)
		{
			foundIndex = root._numOfCurrentStoredObjects - 1 - i + 1;
			break;
		}
	}

	if ((foundIndex != 0) && (root.nodesVals[foundIndex - 1] == val)) // 1, 2
	{
		foundIndex--;

		if (root.isLeaf == true) // 1
		{
			for (uint16_t i = foundIndex; i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				root.nodesVals[i] = root.nodesVals[i + 1];

			root._numOfCurrentStoredObjects--;

			saveNode(root);
		}
		else // 2
		{
			Node<T> leftNeighbourNode(maxNumOfObjectsInNode, root._levelMarker - 1), rightNeighbourNode(maxNumOfObjectsInNode, root._levelMarker - 1);

			getNode(root.childrenNodesIds[foundIndex], leftNeighbourNode);		// left
			getNode(root.childrenNodesIds[foundIndex + 1], rightNeighbourNode); // right

			if (leftNeighbourNode._numOfCurrentStoredObjects >= t) // a
			{
				while (leftNeighbourNode.isLeaf != true)
					getNode(leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects], leftNeighbourNode);

				T valToMove = leftNeighbourNode.nodesVals[leftNeighbourNode._numOfCurrentStoredObjects - 1];

				root.nodesVals[foundIndex] = valToMove;

				saveNode(root);

				getNode(root.childrenNodesIds[foundIndex], leftNeighbourNode);
				deleteKey(leftNeighbourNode, valToMove);
			}
			else if (rightNeighbourNode._numOfCurrentStoredObjects >= t) // b
			{
				while (rightNeighbourNode.isLeaf != true)
					getNode(rightNeighbourNode.childrenNodesIds[0], rightNeighbourNode);

				T valToMove = rightNeighbourNode.nodesVals[0];

				root.nodesVals[foundIndex] = valToMove;

				saveNode(root);

				getNode(root.childrenNodesIds[foundIndex + 1], rightNeighbourNode);
				deleteKey(rightNeighbourNode, valToMove);
			}
			else // c
			{
				leftNeighbourNode.nodesVals[leftNeighbourNode._numOfCurrentStoredObjects] = val;
				leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects + 1] = rightNeighbourNode.childrenNodesIds[0];

				for (size_t i = 0; i < rightNeighbourNode._numOfCurrentStoredObjects; i++)
				{
					leftNeighbourNode.nodesVals[leftNeighbourNode._numOfCurrentStoredObjects + 1 + i] = rightNeighbourNode.nodesVals[i];
					leftNeighbourNode.childrenNodesIds[leftNeighbourNode._numOfCurrentStoredObjects + 1 + 1 + i] = rightNeighbourNode.childrenNodesIds[i + 1];
				}

				leftNeighbourNode._numOfCurrentStoredObjects += (1 + rightNeighbourNode._numOfCurrentStoredObjects);

				for (uint16_t i = foundIndex; i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				{
					root.nodesVals[i] = root.nodesVals[i + 1];
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
		Node<T> mayBeContain(maxNumOfObjectsInNode, root._levelMarker - 1);

		uint16_t mayBeContainRefNumber = (foundIndex != 0) ? (foundIndex - 1 + 1) : 0; // = foundIndex

		getNode(root.childrenNodesIds[mayBeContainRefNumber], mayBeContain);

		// foundIndex--;

		if (mayBeContain._numOfCurrentStoredObjects == minNumOfObjectsInNode)
		{
			Node<T> neigbourOfMayBeContain(maxNumOfObjectsInNode, root._levelMarker - 1);

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
					mayBeContain.nodesVals[mayBeContain._numOfCurrentStoredObjects] = root.nodesVals[mayBeContainRefNumber];
					mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1] = neigbourOfMayBeContain.childrenNodesIds[0];
					mayBeContain._numOfCurrentStoredObjects++;

					root.nodesVals[mayBeContainRefNumber] = neigbourOfMayBeContain.nodesVals[0];

					neigbourOfMayBeContain._numOfCurrentStoredObjects--;

					for (size_t i = 0; i < neigbourOfMayBeContain._numOfCurrentStoredObjects; i++)
					{
						neigbourOfMayBeContain.nodesVals[i] = neigbourOfMayBeContain.nodesVals[i + 1];
						neigbourOfMayBeContain.childrenNodesIds[i] = neigbourOfMayBeContain.childrenNodesIds[i + 1];
					}

					neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects] = neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects + 1];
				}
				else
				{
					for (size_t i = 0; i < mayBeContain._numOfCurrentStoredObjects; i++)
					{
						mayBeContain.nodesVals[mayBeContain._numOfCurrentStoredObjects - 1 - i + 1] = mayBeContain.nodesVals[mayBeContain._numOfCurrentStoredObjects - 1 - i];
						mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1 - 1 - i + 1] = mayBeContain.childrenNodesIds[mayBeContain._numOfCurrentStoredObjects + 1 - 1 - i];
					}
					mayBeContain.childrenNodesIds[1] = mayBeContain.childrenNodesIds[0];

					mayBeContain.nodesVals[0] = root.nodesVals[neigbourOfMayBeContainRefNumber];
					mayBeContain.childrenNodesIds[0] = neigbourOfMayBeContain.childrenNodesIds[neigbourOfMayBeContain._numOfCurrentStoredObjects];
					mayBeContain._numOfCurrentStoredObjects++;

					root.nodesVals[neigbourOfMayBeContainRefNumber] = neigbourOfMayBeContain.nodesVals[neigbourOfMayBeContain._numOfCurrentStoredObjects - 1];

					neigbourOfMayBeContain._numOfCurrentStoredObjects--;
				}

				saveNode(root);
				saveNode(neigbourOfMayBeContain);
				saveNode(mayBeContain);

				deleteKey(mayBeContain, val);
			}
			else // b
			{
				Node<T> *left;
				Node<T> *right;

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

				left->nodesVals[left->_numOfCurrentStoredObjects] = root.nodesVals[std::min(mayBeContainRefNumber, neigbourOfMayBeContainRefNumber)];
				left->childrenNodesIds[left->_numOfCurrentStoredObjects + 1] = right->childrenNodesIds[0];

				for (size_t i = 0; i < right->_numOfCurrentStoredObjects; i++)
				{
					left->nodesVals[left->_numOfCurrentStoredObjects + 1 + i] = right->nodesVals[i];
					left->childrenNodesIds[left->_numOfCurrentStoredObjects + 1 + 1 + i] = right->childrenNodesIds[i + 1];
				}

				left->_numOfCurrentStoredObjects += (1 + right->_numOfCurrentStoredObjects);

				for (uint16_t i = std::min(mayBeContainRefNumber, neigbourOfMayBeContainRefNumber); i <= root._numOfCurrentStoredObjects - 1 - 1; i++)
				{
					root.nodesVals[i] = root.nodesVals[i + 1];
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
							Node<T> newRoot(maxNumOfObjectsInNode, root._levelMarker - 1);
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

template <class T>
B_Tree<T>::~B_Tree()
{
	nodesFile.close();

	if (_keepOnDiskAfterDestruction == false)
		std::filesystem::remove(nodesFileName);
}