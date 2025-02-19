#include "btree.h"

#include <iostream>
#include <vector>
#include <gtest/gtest.h>

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

void createRandomShuffle(uint64_t minVal, uint64_t maxVal, size_t numberOfElements, std::vector<uint64_t> &dst)
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

TEST(BtreeTests, TestNodesSavingGetting)
{
	size_t pageSizeInBytes = 64;

	B_Tree<uint64_t> tr(false, false, std::string("nodes"), pageSizeInBytes, 0);
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

	ASSERT_EQ(nd1.childrenNodesIds.size(), nd2.childrenNodesIds.size());
	ASSERT_EQ(nd1.nodesVals.size(), nd2.nodesVals.size());

	EXPECT_EQ(nd1.isLeaf, nd2.isLeaf);
	EXPECT_EQ(nd1._id, nd2._id);
	EXPECT_EQ(nd1._numOfCurrentStoredObjects, nd2._numOfCurrentStoredObjects);
	EXPECT_EQ(nd1._levelMarker, nd2._levelMarker);

	for (size_t i = 0; i < nd1.nodesVals.size(); i++)
		EXPECT_EQ(nd1.nodesVals[i], nd2.nodesVals[i]);

	for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
		EXPECT_EQ(nd1.childrenNodesIds[i], nd2.childrenNodesIds[i]);
}

TEST(BtreeTests, TestSomeValues)
{
	size_t pageSizeInBytes = 64;
	B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 0);

	std::vector<uint64_t> vals = {9, 10, 3, 134, 104, 17};

	for (auto val : vals)
		btr.insertKey(val);

	for (auto val : vals)
	{
		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		ASSERT_NE(index, 0);
		ASSERT_EQ(findNode.nodesVals[index - 1], val);
	}
}

TEST(BtreeTests, TestRandomValues)
{
	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 100'000, 40, vals);

	size_t pageSizeInBytes = 128;
	B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 4);

	for (auto val : vals)
		btr.insertKey(val);

	for (auto val : vals)
	{
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, val, findNode, index);

		ASSERT_NE(index, 0);
		ASSERT_EQ(findNode.nodesVals[index - 1], val);
	}
}

TEST(BtreeTests, TestMoreValues)
{
	srand(std::time(NULL));
	std::vector<uint64_t> sh;

	createRandomShuffle(1, 100'000, 2000, sh);

	size_t pageSizeInBytes = 4096;
	B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 4);

	for (auto i : sh)
		btr.insertKey(i);

	for (auto i : sh)
	{
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, i, findNode, index);

		ASSERT_NE(index, 0);
		ASSERT_EQ(findNode.nodesVals[index - 1], i);
	}
}

TEST(BtreeTests, TestDeletionOfSomeValues)
{
	size_t pageSizeInBytes = 64;
	B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 0);

	std::vector<uint64_t> vals = {9, 10, 3, 134, 104, 17};

	for (auto val : vals)
		btr.insertKey(val);

	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		ASSERT_EQ(index, 0);
	}
}

TEST(BtreeTests, TestDeletionOfRandomValues)
{
	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 100'000, 40, vals);

	size_t pageSizeInBytes = 128;
	B_Tree<uint64_t> btr(false, true, std::string("nodes"), pageSizeInBytes, 4);

	for (auto val : vals)
		btr.insertKey(val);

	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		ASSERT_EQ(index, 0);
	}
}

TEST(BtreeTests, TestDeletionOfMoreValues)
{
	srand(std::time(NULL));
	std::vector<uint64_t> vals;

	createRandomShuffle(1, 1'000'000, 15'000, vals);

	size_t pageSizeInBytes = 4096;
	B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 5);

	for (auto val : vals)
		btr.insertKey(val);

	for (auto val : vals)
	{
		btr.deleteKey(btr._root, val);

		uint64_t searchVal = val;
		Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

		uint16_t index = 0;
		btr.search(btr._root, searchVal, findNode, index);

		ASSERT_EQ(index, 0);
	}
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}