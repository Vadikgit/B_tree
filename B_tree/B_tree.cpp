// B_tree.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _CRT_SECURE_NO_WARNINGS
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

void createRandomShuffle(uint64_t minVal, uint64_t maxVal, size_t numberOfElements, std::vector<uint64_t>& dst)
{	
	dst.resize(maxVal - minVal + 1);

	for (uint64_t i = minVal; i <= maxVal; i++)
		dst[i - minVal] = i;

	for (size_t i = 0; i < dst.size() - 1; i++)
	{
		size_t pos = rand() % (dst.size() - i - 1) + i + 1;
		std::swap(dst[i], dst[pos]);
	}

	dst.resize(numberOfElements);
}



class uint32char64type
{
public:
	uint32_t _key;
	char _val[128];
};

bool operator<(const uint32char64type& tp1, const uint32char64type& tp2)
{
	return (tp1._key < tp2._key);
}
bool operator==(const uint32char64type& tp1, const uint32char64type& tp2)
{
	return tp1._key == tp2._key;
}
bool operator<=(const uint32char64type& tp1, const uint32char64type& tp2)
{
	return tp1._key <= tp2._key;
}
bool operator>(const uint32char64type& tp1, const uint32char64type& tp2)
{
	return tp1._key > tp2._key;
}

template <class T>
class Node
{
public:
	Node(uint16_t numOfObjectsInNode, uint32_t levelMarker) : _levelMarker{ levelMarker }, _id{ 0 }, isLeaf{ false }, _numOfCurrentStoredObjects{ 0 }
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

template <class T>
class B_Tree;

struct levelMarkersTreeNode
{
	uint32_t levelMarker;
	uint32_t id;
};

auto levelMarkersTreeNodeComp = [](const levelMarkersTreeNode& n1, const levelMarkersTreeNode& n2) -> bool {
		return ((n1.levelMarker < n2.levelMarker) || (n1.levelMarker == n2.levelMarker && n1.id < n2.id));
};

template <class T>
class NodeCache
{
public:
	size_t _maxNumOfPages;
	std::unordered_map<uint32_t, Node<T>> _buffer;
	std::set <levelMarkersTreeNode,  decltype(levelMarkersTreeNodeComp)> _levelMarkersTree;
	B_Tree<T>* _ptrToTree;

	NodeCache(B_Tree<T>* ptrToTree, size_t maxNumOfPages) : _ptrToTree{ ptrToTree }, _maxNumOfPages { maxNumOfPages }, _levelMarkersTree(levelMarkersTreeNodeComp) {};

	bool processPretendentToCache(Node<T>& node, bool getSave);
	void getNode(uint32_t id, Node<T>& dst);
	void saveNode(Node<T>& src);
	void freeNode(Node<T>& node);

	void writeNodeToDisk(Node<T>& src);
	void readNodeFromDisk(uint32_t id, Node<T>& dst);
};

template <class T>
bool NodeCache<T>::processPretendentToCache(Node<T>& node, bool getSave)
{
	if (_buffer.size() != _maxNumOfPages)
	{
		_buffer.insert({ node._id, node });
		_levelMarkersTree.insert({ node._levelMarker, node._id });

		return true;
	}

	auto minimum = *(_levelMarkersTree.begin());

	if (getSave == true)
	{
		if (minimum.levelMarker <= node._levelMarker)
		{
			_levelMarkersTree.erase(minimum);
			_levelMarkersTree.insert({ node._levelMarker, node._id });

			writeNodeToDisk(_buffer.find(minimum.id)->second);

			_buffer.erase(minimum.id);
			_buffer.insert({ node._id, node });

			return true;
		}
	}
	else
	{
		if (minimum.levelMarker <= node._levelMarker)
		{
			_levelMarkersTree.erase(minimum);
			_levelMarkersTree.insert({ node._levelMarker, node._id });

			writeNodeToDisk(_buffer.find(minimum.id)->second);

			_buffer.erase(minimum.id);
			_buffer.insert({ node._id, node });

			return true;
		}
	}

	return false;
}

template <class T>
void NodeCache<T>::readNodeFromDisk(uint32_t id, Node<T>& dst)
{

	/*std::chrono::time_point<std::chrono::system_clock> t1, t2;

	dst._id = id; // id

	std::stringstream stream;
	stream << std::hex << id;
	std::string nodeFileName(stream.str());
	nodeFileName += ".btrnd";
	nodeFileName = _ptrToTree->pathToNodes + std::string("\\") + nodeFileName;

	std::vector<uint8_t> nodeBytes;
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	//std::ifstream nodeFile;
	//nodeFile.open(nodeFileName, std::ios::binary | std::ios::in);
	//nodeFile.read((char*)&(nodeBytes[0]), _ptrToTree->_fileSizeInBytes);

	//t1 = std::chrono::system_clock::now();
	FILE* input_file = fopen(nodeFileName.c_str(), "rb");
	//t2 = std::chrono::system_clock::now();

	t1 = std::chrono::system_clock::now();
	fread((char*)&(nodeBytes[0]), nodeBytes.size(), 1, input_file);
	t2 = std::chrono::system_clock::now();

	fclose(input_file);
	//t2 = std::chrono::system_clock::now();
	////////////////////////////////

	size_t numOfPassedBytes = 0;

	dst.isLeaf = *((bool*)&(nodeBytes[numOfPassedBytes])); // isLeaf
	numOfPassedBytes += sizeof(bool);

	dst._numOfCurrentStoredObjects = *((uint16_t*)&(nodeBytes[numOfPassedBytes])); // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	dst._levelMarker = *((uint32_t*)&(nodeBytes[numOfPassedBytes])); // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	//for (size_t i = 0; i < dst.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	dst.childrenNodesIds[i] = *((uint32_t*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(uint32_t);
	//}
	std::memcpy((void*)&(dst.childrenNodesIds[0]), (void*)&(nodeBytes[numOfPassedBytes]), dst.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * dst.childrenNodesIds.size();

	//for (size_t i = 0; i < dst.nodesVals.size(); i++) // nodesVals
	//{
	//	dst.nodesVals[i] = *((T*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(T);
	//}
	std::memcpy((void*)&(dst.nodesVals[0]), (void*)&(nodeBytes[numOfPassedBytes]), dst.nodesVals.size() * sizeof(dst.nodesVals[0]));
	numOfPassedBytes += sizeof(dst.nodesVals[0]) * dst.nodesVals.size();

	

	_ptrToTree->_rdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_rnum++;*/



	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	dst._id = id; // id

	std::vector<uint8_t> nodeBytes;
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	t1 = std::chrono::system_clock::now();
	_ptrToTree->nodesFile.seekg(id * _ptrToTree->_fileSizeInBytes);
	_ptrToTree->nodesFile.read((char*)&(nodeBytes[0]), nodeBytes.size());
	t2 = std::chrono::system_clock::now();

	
	////////////////////////////////

	size_t numOfPassedBytes = 0;

	dst.isLeaf = *((bool*)&(nodeBytes[numOfPassedBytes])); // isLeaf
	numOfPassedBytes += sizeof(bool);

	dst._numOfCurrentStoredObjects = *((uint16_t*)&(nodeBytes[numOfPassedBytes])); // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	dst._levelMarker = *((uint32_t*)&(nodeBytes[numOfPassedBytes])); // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	//for (size_t i = 0; i < dst.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	dst.childrenNodesIds[i] = *((uint32_t*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(uint32_t);
	//}
	std::memcpy((void*)&(dst.childrenNodesIds[0]), (void*)&(nodeBytes[numOfPassedBytes]), dst.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * dst.childrenNodesIds.size();

	//for (size_t i = 0; i < dst.nodesVals.size(); i++) // nodesVals
	//{
	//	dst.nodesVals[i] = *((T*)&(nodeBytes[numOfPassedBytes]));
	//	numOfPassedBytes += sizeof(T);
	//}
	std::memcpy((void*)&(dst.nodesVals[0]), (void*)&(nodeBytes[numOfPassedBytes]), dst.nodesVals.size() * sizeof(dst.nodesVals[0]));
	numOfPassedBytes += sizeof(dst.nodesVals[0]) * dst.nodesVals.size();



	_ptrToTree->_rdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_rnum++;
}

template <class T>
void NodeCache<T>::writeNodeToDisk(Node<T>& src)
{
	/*std::chrono::time_point<std::chrono::system_clock> t1, t2;

	size_t numOfPassedBytes = 0;
	std::vector<uint8_t> nodeBytes;
//	nodeBytes.assign((sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + sizeof(src.nodesVals[0]) * _ptrToTree->maxNumOfObjectsInNode + sizeof(src.childrenNodesIds[0]) * (_ptrToTree->maxNumOfObjectsInNode + 1)), 0);
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	*((bool*)&(nodeBytes[numOfPassedBytes])) = src.isLeaf; // isLeaf
	numOfPassedBytes += sizeof(bool);

	*((uint16_t*)&(nodeBytes[numOfPassedBytes])) = src._numOfCurrentStoredObjects; // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	*((uint32_t*)&(nodeBytes[numOfPassedBytes])) = src._levelMarker; // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	//for (size_t i = 0; i < src.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	*((uint32_t*)&(nodeBytes[numOfPassedBytes])) = src.childrenNodesIds[i];
	//	numOfPassedBytes += sizeof(uint32_t);
	//}
	std::memcpy((void*)&(nodeBytes[numOfPassedBytes]), (void*)&(src.childrenNodesIds[0]), src.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * src.childrenNodesIds.size();

	//for (size_t i = 0; i < src.nodesVals.size(); i++) // nodesVals
	//{
	//	*((T*)&(nodeBytes[numOfPassedBytes])) = src.nodesVals[i];
	//	numOfPassedBytes += sizeof(T);
	//}
	std::memcpy((void *)&(nodeBytes[numOfPassedBytes]), (void *)&(src.nodesVals[0]), src.nodesVals.size() * sizeof(src.nodesVals[0]));
	numOfPassedBytes += sizeof(T) * src.nodesVals.size();

	////////////////////////////////

	std::stringstream stream;
	stream << std::hex << src._id;
	std::string nodeFileName(stream.str());
	nodeFileName += ".btrnd";
	nodeFileName = _ptrToTree->pathToNodes + std::string("\\") + nodeFileName;

	//std::ofstream nodeFile;
	//nodeFile.open(nodeFileName, std::ios::binary | std::ios::trunc | std::ios::out);
	//nodeFile.write((char*)&(nodeBytes[0]), nodeBytes.size());
	
	//t1 = std::chrono::system_clock::now();
	FILE* output_file = fopen(nodeFileName.c_str(), "wb");
	//t2 = std::chrono::system_clock::now();

	t1 = std::chrono::system_clock::now();
	fwrite((char*)&(nodeBytes[0]), nodeBytes.size(), 1, output_file);
	t2 = std::chrono::system_clock::now();

	//fflush(output_file);
	//t2 = std::chrono::system_clock::now();
	fclose(output_file);
	//t2 = std::chrono::system_clock::now();

	_ptrToTree->_wdur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	_ptrToTree->_wnum++;*/


	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	size_t numOfPassedBytes = 0;
	std::vector<uint8_t> nodeBytes;
	//	nodeBytes.assign((sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + sizeof(src.nodesVals[0]) * _ptrToTree->maxNumOfObjectsInNode + sizeof(src.childrenNodesIds[0]) * (_ptrToTree->maxNumOfObjectsInNode + 1)), 0);
	nodeBytes.assign(_ptrToTree->_fileSizeInBytes, 0);

	*((bool*)&(nodeBytes[numOfPassedBytes])) = src.isLeaf; // isLeaf
	numOfPassedBytes += sizeof(bool);

	*((uint16_t*)&(nodeBytes[numOfPassedBytes])) = src._numOfCurrentStoredObjects; // numOfCurrentStoredObjects
	numOfPassedBytes += sizeof(uint16_t);

	*((uint32_t*)&(nodeBytes[numOfPassedBytes])) = src._levelMarker; // levelMarker
	numOfPassedBytes += sizeof(uint32_t);

	//for (size_t i = 0; i < src.childrenNodesIds.size(); i++) // childrenNodesIds
	//{
	//	*((uint32_t*)&(nodeBytes[numOfPassedBytes])) = src.childrenNodesIds[i];
	//	numOfPassedBytes += sizeof(uint32_t);
	//}
	std::memcpy((void*)&(nodeBytes[numOfPassedBytes]), (void*)&(src.childrenNodesIds[0]), src.childrenNodesIds.size() * sizeof(uint32_t));
	numOfPassedBytes += sizeof(uint32_t) * src.childrenNodesIds.size();

	//for (size_t i = 0; i < src.nodesVals.size(); i++) // nodesVals
	//{
	//	*((T*)&(nodeBytes[numOfPassedBytes])) = src.nodesVals[i];
	//	numOfPassedBytes += sizeof(T);
	//}
	std::memcpy((void*)&(nodeBytes[numOfPassedBytes]), (void*)&(src.nodesVals[0]), src.nodesVals.size() * sizeof(src.nodesVals[0]));
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
void NodeCache<T>::getNode(uint32_t id, Node<T>& dst)
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
void NodeCache<T>::saveNode(Node<T>& src)
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
void NodeCache<T>::freeNode(Node<T>& node)
{
	// remove from _freqTree and from buffer
	if (_buffer.find(node._id) != _buffer.end())
	{
		_levelMarkersTree.erase({ node._levelMarker, node._id });
		_buffer.erase(node._id);
	}

	// delete from disk
	/*std::stringstream stream;
	stream << std::hex << node._id;
	std::string nodeFileName(stream.str());
	nodeFileName += ".btrnd";
	nodeFileName = _ptrToTree->pathToNodes + std::string("\\") + nodeFileName;

	std::filesystem::remove(nodeFileName);*/


}



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
	uint16_t minNumOfObjectsInNode; // = t - 1
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

	void getNode(uint32_t id, Node<T>& dst);
	void saveNode(Node<T>& src);
	void initNode(Node<T>& node)
	{
		node._id = largestExistingId + 1;
		largestExistingId += 1;
	}
	void freeNode(Node<T>& node);
	void makeRootUpper();
	void makeRootLower(Node<T>& newRoot);
	void createTree();
	void search(Node<T>& root, T& val, Node<T>& dst_node, uint16_t& indexWhereValueFound);
	void splitChild(Node<T>& node, uint16_t index);
	void insertKey(T& val);
	void insertIntoNonFull(Node<T>& node, T& val);
	void deleteKey(Node<T>& root, T& val);
};

template <class T>
B_Tree<T>::B_Tree(bool recover, bool keepOnDiskAfterDestruction, std::string _pathToNodes, size_t diskPageSizeInBytes, size_t maxNumOfPagesInCache) :
	maxNumOfObjectsInNode{ uint16_t((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) 
	- ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) + 1) % 2 )},
	_root((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t))
	- ((diskPageSizeInBytes - sizeof(uint32_t) - sizeof(uint16_t) - sizeof(bool) - sizeof(uint32_t)) / (sizeof(T) + sizeof(uint32_t)) + 1) % 2, 0),
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
	_fileSizeInBytes = sizeof(bool) + sizeof(uint16_t) + sizeof(uint32_t) + sizeof(T) * maxNumOfObjectsInNode 
					 + sizeof(uint32_t) * (maxNumOfObjectsInNode + 1);
	
	//////
	nodesFileName = pathToNodes + std::string("\\file.btr");
	std::ofstream(nodesFileName.c_str()).close();
	nodesFile.open(nodesFileName.c_str(), std::ios::binary | std::ios::out | std::ios::in);
	if (!nodesFile)
	{
		std::cerr << "\"file.btr\" could not be opened!" << std::endl;
		exit(1);
	}
	//////

	std::cout << "\n<< " << maxNumOfObjectsInNode << " objects in node>>\n";
	createTree();
}

template <class T>
void B_Tree<T>::getNode(uint32_t id, Node<T>& dst)
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
void B_Tree<T>::saveNode(Node<T>& src)
{
	if (src._id == _root._id)
		_root = src;
	else
	{
		_cache.saveNode(src);
	}
}

template <class T>
void B_Tree<T>::freeNode(Node<T>& node)
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
void B_Tree<T>::makeRootLower(Node<T>& newRoot)
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
void B_Tree<T>::search(Node<T>& root, T& val, Node<T>& dst_node, uint16_t& indexWhereValueFound)
{
	uint16_t i = 1;

	//while ((i <= root._numOfCurrentStoredObjects) && (val > root.nodesVals[i - 1]))
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
void B_Tree<T>::splitChild(Node<T>& node, uint16_t index)
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
void B_Tree<T>::insertIntoNonFull(Node<T>& node, T& val) 
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
void B_Tree<T>::insertKey(T& val)
{
	if (_root._numOfCurrentStoredObjects == maxNumOfObjectsInNode)
	{
		makeRootUpper();
		splitChild(_root, 1); // previos root will be saved here
	}

	insertIntoNonFull(_root, val);
}

template <class T>
void B_Tree<T>::deleteKey(Node<T>& root, T& val)
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

			getNode(root.childrenNodesIds[foundIndex], leftNeighbourNode); // left
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

		//foundIndex--;

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
				Node<T>* left;
				Node<T>* right;

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
/*	if (_keepOnDiskAfterDestruction == false)
	{
		for (const auto& entry : std::filesystem::directory_iterator(pathToNodes))
			std::filesystem::remove(entry.path());
	}*/

	nodesFile.close();

	if (_keepOnDiskAfterDestruction == false)
		std::filesystem::remove(nodesFileName);
}


//-----------------------------------------------------
void testNodesSavingGetting(std::string path)
{
	std::cout << "\n\n=============== testNodesSavingGetting ==============\n";

	size_t pageSizeInBytes = 64;

	B_Tree<uint64_t> tr(false, false, path, pageSizeInBytes, 0);
	Node<uint64_t> nd1(tr.maxNumOfObjectsInNode, 45);
	tr.initNode(nd1);

	nd1._numOfCurrentStoredObjects = nd1.nodesVals.size();

	for (size_t i = 0; i < nd1.nodesVals.size(); i++)
		nd1.nodesVals[i] = (i * 9487 + 234);

	for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
		nd1.childrenNodesIds[i] = (i * 9875 + 34);

	tr._cache.saveNode(nd1);

	Node<uint64_t> nd2(tr.maxNumOfObjectsInNode, 0);
	tr._cache.getNode(2, nd2);

	if (nd2.childrenNodesIds.size() != nd1.childrenNodesIds.size())
	{
		std::cout << "children nodes ids arrs sizes not equal";
		return ;
	}

	if (nd2.nodesVals.size() != nd1.nodesVals.size())
	{
		std::cout << "nodes val arrs sizes not equal";
		return ;
	}


	std::cout << "isLeaf: " << nd1.isLeaf << ' ' << nd2.isLeaf << "\n";
	std::cout << "_id: " << nd1._id << ' ' << nd2._id << "\n";
	std::cout << "_numOfCurrentStoredObjects: " << nd1._numOfCurrentStoredObjects << ' ' << nd2._numOfCurrentStoredObjects << "\n";
	std::cout << "_levelMarker: " << nd1._levelMarker << ' ' << nd2._levelMarker << "\n";


	std::cout << "\n\n";

	std::cout << "vals: \n";
	for (size_t i = 0; i < nd1.nodesVals.size(); i++)
		std::cout << nd1.nodesVals[i] << ' ' << nd2.nodesVals[i] << '\n';

	std::cout << "\n\n";

	std::cout << "childs: \n";
	for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
		std::cout << nd1.childrenNodesIds[i] << ' ' << nd2.childrenNodesIds[i] << '\n';

	std::cout << "\n==============================================================\n\n";
}

void testUserTypeNodesSavingGetting(std::string path)
{
	std::cout << "\n\n=============== testUserTypeNodesSavingGetting ==============\n";

	size_t pageSizeInBytes = 1024;

	B_Tree<uint32char64type> tr(false, false, path, pageSizeInBytes, 0);
	Node<uint32char64type> nd1(tr.maxNumOfObjectsInNode, 45);
	tr.initNode(nd1);

	nd1._numOfCurrentStoredObjects = nd1.nodesVals.size();

	for (size_t i = 0; i < nd1.nodesVals.size(); i++)
	{
		nd1.nodesVals[i]._key = (i * 9487 + 234);

		for (size_t j = 0; j < sizeof(nd1.nodesVals[i]._val); j++)
			nd1.nodesVals[i]._val[j] = i + j;
	}
		
	for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
		nd1.childrenNodesIds[i] = (i * 9875 + 34);

	tr._cache.saveNode(nd1);

	Node<uint32char64type> nd2(tr.maxNumOfObjectsInNode, 0);
	tr._cache.getNode(2, nd2);

	if (nd2.childrenNodesIds.size() != nd1.childrenNodesIds.size())
	{
		std::cout << "children nodes ids arrs sizes not equal";
		return;
	}

	if (nd2.nodesVals.size() != nd1.nodesVals.size())
	{
		std::cout << "nodes val arrs sizes not equal";
		return;
	}


	std::cout << "isLeaf: " << nd1.isLeaf << ' ' << nd2.isLeaf << "\n";
	std::cout << "_id: " << nd1._id << ' ' << nd2._id << "\n";
	std::cout << "_numOfCurrentStoredObjects: " << nd1._numOfCurrentStoredObjects << ' ' << nd2._numOfCurrentStoredObjects << "\n";
	std::cout << "_levelMarker: " << nd1._levelMarker << ' ' << nd2._levelMarker << "\n";


	std::cout << "\n\n";

	std::cout << "vals: \n";
	for (size_t i = 0; i < nd1.nodesVals.size(); i++)
	{
		std::cout << "\n\t" << i << "\n";
		for (size_t j = 0; j < sizeof(nd1.nodesVals[i]._val); j++)
		{
			std::cout << (int) nd1.nodesVals[i]._val[j] << ' ' << (int)nd2.nodesVals[i]._val[j] << '\n';
		}
		std::cout << "\n";
	}

	std::cout << "\n\n";

	std::cout << "childs: \n";
	for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
		std::cout << nd1.childrenNodesIds[i] << ' ' << nd2.childrenNodesIds[i] << '\n';

	std::cout << "\n==============================================================\n\n";
}

void testSomeValues(std::string path)
{
	std::cout << "\n\n=============== testSomeValues ==============\n";

	size_t pageSizeInBytes = 64;
	B_Tree<uint64_t> btr(false, false, path, pageSizeInBytes, 0);

	std::vector<uint64_t> vals = { 9, 10, 3, 134, 104, 17 };

	for (auto val : vals)
		btr.insertKey(val);
	
	for (auto val : vals)
	{
		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		if (index == 0)
		{
			std::cout << "not founded(\n";
		}
		else
		{
			std::cout << "index: " << index << '\n';
			for (size_t i = 0; i < findNode.nodesVals.size(); i++)
				std::cout << findNode.nodesVals[i] << ' ';
			std::cout << '\n';
		}
	}

	std::cout << "\n==============================================================\n\n";
}

void testRandomValues(std::string path)
{
	std::cout << "\n\n=============== testRandomValues ==============\n";

	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 100'000, 40, vals);


	size_t pageSizeInBytes = 128;
	B_Tree<uint64_t> btr(false, false, "nodes", pageSizeInBytes, 4);

	for (auto val : vals)
		btr.insertKey(val);

	std::cout << "filled\n";

	for (auto val : vals) {

		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, val, findNode, index);

		if (index == 0)
			std::cout << "\n\nvalue " << val << " not found(\n";
		else
		{
			std::cout << "\n\nvalue " << val << " found\n";
			std::cout << "index: " << index << '\n';
			for (size_t i = 0; i < findNode.nodesVals.size(); i++)
			{
				std::cout << findNode.nodesVals[i] << ' ';
			}
		}
	}

	std::cout << "\n==============================================================\n\n";
}

void testMoreValues(std::string path)
{
	std::cout << "\n\n=============== testMoreValues ==============\n";
	srand(std::time(NULL));
	std::vector<uint64_t> sh;

	createRandomShuffle(1, 100'000, 2000, sh);


	size_t pageSizeInBytes = 4096;
	B_Tree<uint64_t> btr(false, false, "nodes", pageSizeInBytes, 4);

	for (auto i : sh)
		btr.insertKey(i);

	std::cout << "filled\n";

	for (auto i : sh) {

		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, i, findNode, index);

		if (index == 0)
			std::cout << "value " << i << " not found(\n";

	}

	std::cout << "\n==============================================================\n\n";
}

void testDeletionOfSomeValues(std::string path)
{
	std::cout << "\n\n=============== testDeletionOfSomeValues ==============\n";

	size_t pageSizeInBytes = 64;
	B_Tree<uint64_t> btr(false, false, path, pageSizeInBytes, 0);

	std::vector<uint64_t> vals = { 9, 10, 3, 134, 104, 17 };

	for (auto val : vals)
		btr.insertKey(val);

	std::cout << "\nfilled\n";


	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		std::cout << "\n" << val << " deleted\n";

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		if (index == 0)
		{
			std::cout << "not founded(\n";
		}
		else
		{
			std::cout << "index: " << index << '\n';
			for (size_t i = 0; i < findNode.nodesVals.size(); i++)
				std::cout << findNode.nodesVals[i] << ' ';
			std::cout << '\n';
		}
	}

	std::cout << "\n==============================================================\n\n";
}

void testDeletionOfRandomValues(std::string path)
{
	std::cout << "\n\n=============== testDeletionOfRandomValues ==============\n";

	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 100'000, 40, vals);


	size_t pageSizeInBytes = 128;
	B_Tree<uint64_t> btr(false, true, "nodes", pageSizeInBytes, 4);

	for (auto val : vals)
		btr.insertKey(val);

	std::cout << "\nfilled\n";

	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		std::cout << "\n" << val << " deleted\n";

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		if (index == 0)
		{
			std::cout << "not founded(\n";
		}
		else
		{
			std::cout << "index: " << index << '\n';
			for (size_t i = 0; i < findNode.nodesVals.size(); i++)
				std::cout << findNode.nodesVals[i] << ' ';
			std::cout << '\n';
		}
	}

	std::cout << "\n==============================================================\n\n";
}

void testDeletionOfMoreValues(std::string path)
{
	std::cout << "\n\n=============== testDeletionOfMoreValues ==============\n";

	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 1'000'000, 15'000, vals);


	size_t pageSizeInBytes = 4096;
	B_Tree<uint64_t> btr(false, false, "nodes", pageSizeInBytes, 5);

	for (auto val : vals)
		btr.insertKey(val);

	std::cout << "\nfilled\n";

	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		//std::cout << "\n" << val << " deleted\n";

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		if (index != 0)
			std::cout << "value " << val << " was found, but deleted(\n";
	}

	std::cout << "\n==============================================================\n\n";
}

void testTimeOfInsertionDeletionOfMoreValues(std::string path)
{
	std::cout << "\n\n=============== testTimeOfInsertionDeletionOfMoreValues ==============\n";

	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	size_t pageSizeInBytes = 4096;
	B_Tree<uint64_t> btr(false, false, "nodes", pageSizeInBytes, 300);


	std::chrono::microseconds dur(0);
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	size_t numOfValuesToProcess = 1'000'000;
	size_t numOfValuesInOneSubShuffle = 2000;
	//size_t numOfPossibleValuesForOneSubhuffle = 100'000;

	std::vector<uint64_t> partitionIndexesShuffle;
	createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // insertion
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
			btr.insertKey(val);

		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nfilled\n" << dur.count() << " mcs;\tread: " << btr._rdur.count() << " mcs, " << btr._rnum << " times;\twrite: " << btr._wdur.count() << " mcs, " << btr._wnum << " times" << std::endl;


	dur = std::chrono::microseconds(0);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
		{
			uint64_t searchVal = val;
			Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

			uint16_t index = 0;
			btr.search(btr._root, searchVal, findNode, index);

			if (index == 0)
				std::cout << "value " << val << " not found(\n";
		}
		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nsearched\n" << dur.count() << " mcs" << std::endl;


	dur = std::chrono::microseconds(0);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search + del
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
		{
			uint64_t searchVal = val;
			Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

			uint16_t index = 0;
			btr.search(btr._root, searchVal, findNode, index);

			if (index == 0)
				std::cout << "value " << val << " not found(\n";

			btr.deleteKey(btr._root, val);
		}
		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nsearched + deleted\n" << dur.count() << " mcs" << std::endl;

	std::cout << "\n==============================================================\n\n";

}

void testTimeOfInsertionDeletionOfMoreUserTypeValues(std::string path)
{
	std::cout << "\n\n=============== testTimeOfInsertionDeletionOfMoreUserTypeValues ==============\n";

	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	size_t pageSizeInBytes = 4096;
	B_Tree<uint32char64type> btr(false, true, "nodes", pageSizeInBytes, 100);


	std::chrono::microseconds dur(0);
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	size_t numOfValuesToProcess = 1'000'000;
	size_t numOfValuesInOneSubShuffle = 2000;
	//size_t numOfPossibleValuesForOneSubhuffle = 100'000;

	std::vector<uint64_t> partitionIndexesShuffle;
	createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // insertion
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
		{
			uint32char64type t;
			t._key = val;
			for (size_t i = 0; i < sizeof(t._val); i++)
				t._val[i] = (val % 10) + i;

			btr.insertKey(t);
		}

		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nfilled\n" << dur.count() << " mcs;\tread: " << btr._rdur.count() << " mcs, " << btr._rnum << " times;\twrite: " << btr._wdur.count() << " mcs, " << btr._wnum << " times" << std::endl;


	dur = std::chrono::microseconds(0);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
		{
			uint32char64type searchVal;
			searchVal._key = val;

			Node<uint32char64type> findNode(btr.maxNumOfObjectsInNode, 0);

			uint16_t index = 0;
			btr.search(btr._root, searchVal, findNode, index);

			if (index == 0)
				std::cout << "value with key" << val << " not found(\n";
		}
		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nsearched\n" << dur.count() << " mcs" << std::endl;


	dur = std::chrono::microseconds(0);

	for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search + del
	{
		createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

		t1 = std::chrono::system_clock::now();

		for (auto val : vals)
		{
			uint32char64type searchVal;
			searchVal._key = val;

			Node<uint32char64type> findNode(btr.maxNumOfObjectsInNode, 0);

			uint16_t index = 0;
			btr.search(btr._root, searchVal, findNode, index);

			if (index == 0)
				std::cout << "value " << val << " not found(\n";

			btr.deleteKey(btr._root, searchVal);
		}
		t2 = std::chrono::system_clock::now();

		dur += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	}
	std::cout << "\nsearched + deleted\n" << dur.count() << " mcs" << std::endl;

	std::cout << "\n==============================================================\n\n";

}

//-------------------------------------------------------


int main()
{
	std::string path = "nodes";

	//testNodesSavingGetting(path);
	//testUserTypeNodesSavingGetting(path);
	//testSomeValues(path);
	//testDeletionOfRandomValues(path);
	//testDeletionOfSomeValues(path);
	//testRandomValues(path);
	//testMoreValues(path);
	//testDeletionOfMoreValues(path);
	testTimeOfInsertionDeletionOfMoreValues(path);
	//testTimeOfInsertionDeletionOfMoreUserTypeValues(path);

	/*uint32char64type t1, t2;
	t1._key = 34;
	for (size_t i = 0; i < sizeof(t1._val); i++)
	{
		t1._val[i] = 97 + i;
	}
	
	t2 = t1;
	t1._val[7] = 97;


	std::cout << sizeof(t1) << ' ' << sizeof(uint32char64type) << ' ' <<(size_t) &(t1._val[0]) <<'\n';
	for (size_t i = 0; i < sizeof(t1._val); i++)
	{
		std::cout << t1._val[i] << ' ';
	}


	std::cout << "\n\n" << sizeof(t2) << ' ' << sizeof(uint32char64type) << ' ' << (size_t) &(t2._val[0]) << '\n';
	for (size_t i = 0; i < sizeof(t2._val); i++)
	{
		std::cout << t2._val[i] << ' ';
	}

	t2._key = 33;

	std::cout << "\n\nt2 less then t1 " << (t2 < t1);*/

	/*std::ofstream("file.bin").close();

	std::fstream fl;
	fl.open("file.bin", std::ios::binary | std::ios::out | std::ios::in);
	if (!fl)
	{
		std::cerr << "Uh oh, SomeText.txt could not be opened for writing!" << std::endl;
		exit(1);
	}


	std::string str = "errhnvu3y7yub8q9nwcvrev657";
	fl.write((char*)&(str[0]), str.length());
	
	std::getchar();

	fl.seekp(100);
	int a = 256 * 256 - 1;
	fl.write((char*)&(a), 2);*/

	/*std::chrono::microseconds dur(0);
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	std::vector<uint64_t> vec;


	vec.resize(100'000);

	for (size_t i = 0; i < vec.size(); i++)
	{
		vec[i] = i;
	}

	t1 = std::chrono::system_clock::now();
	for (size_t i = 0; i < vec.size(); i++)
	{
		for (size_t j = 0; j < vec.size(); j++)
		{
			if (vec[j] == i)
				break;
		}
	}
	t2 = std::chrono::system_clock::now();

	dur = (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	std::cout << "\nlinear search\n" << dur.count() << " mcs" << std::endl;

	t1 = std::chrono::system_clock::now();
	for (size_t i = 0; i < vec.size(); i++)
	{
		std::upper_bound(vec.begin(), vec.end(), i);
	}
	t2 = std::chrono::system_clock::now();

	dur = (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	std::cout << "\nbinary search\n" << dur.count() << " mcs" << std::endl;*/

	/*std::chrono::microseconds dur(0);
	std::chrono::time_point<std::chrono::system_clock> t1, t2;

	std::vector<std::vector<uint64_t>> vec;
	vec.resize(100);

	for (size_t i = 0; i < vec.size(); i++)
	{
		vec[i].resize(100'000);

		for (size_t j = 0; j < vec.size(); j++)
		{
			vec[i][j] = j;
		}
	}

	t1 = std::chrono::system_clock::now();
	auto vec1 = vec;
	
	t2 = std::chrono::system_clock::now();

	dur = (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	std::cout << "\ncopying\n" << dur.count() << " mcs" << std::endl;

	t1 = std::chrono::system_clock::now();
	auto vec2 = std::move(vec1);
	
	t2 = std::chrono::system_clock::now();

	dur = (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
	std::cout << "\nmoving\n" << dur.count() << " mcs" << std::endl;

	std::cout << vec1.empty();*/
}

