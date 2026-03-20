#include "../btree.h"

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

struct uint64_type : public BTProcessable
{
    uint64_t val{0};

    virtual void fromBytes(const std::vector<uint8_t> &ibuf, size_t firstBytePos) override;
    virtual void toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos) override;
    virtual bool less(const BTProcessable &w2) const override;
    virtual size_t sizeInBytes() override;
    virtual uint64_type *createNew() override;

    uint64_type() : BTProcessable() {};
    uint64_type(uint64_t v) : val{v}, BTProcessable() {};
    uint64_type(const uint64_type &other) : val{other.val}, BTProcessable(other) {};
    uint64_type(uint64_type &&other) : val{other.val}, BTProcessable(other) {};

    virtual uint64_type &operator=(const uint64_type &other);
    virtual uint64_type &operator=(uint64_type &&other);

    virtual uint64_type &operator=(const BTProcessable &other) override;
    virtual uint64_type &operator=(BTProcessable &&other) override;
    virtual ~uint64_type() {};
};

void uint64_type::fromBytes(const std::vector<uint8_t> &ibuf, size_t firstBytePos)
{
    val = *(reinterpret_cast<const uint64_t *>(&ibuf[firstBytePos]));
}

void uint64_type::toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos)
{
    *(reinterpret_cast<uint64_t *>(&obuf[firstBytePos])) = val;
}

bool uint64_type::less(const BTProcessable &p2) const
{
    return val < (dynamic_cast<const uint64_type &>(p2)).val;
}

size_t uint64_type::sizeInBytes()
{
    return sizeof(val);
}

uint64_type *uint64_type::createNew()
{
    return new uint64_type();
}

uint64_type &uint64_type::operator=(const uint64_type &other)
{
    if (this != &other)
    {
        val = other.val;
    }

    return *this;
}

uint64_type &uint64_type::operator=(uint64_type &&other)
{
    val = other.val;

    return *this;
}

uint64_type &uint64_type::operator=(const BTProcessable &other)
{
    if (this != &other)
    {
        val = reinterpret_cast<const uint64_type &>(other).val;
    }

    return *this;
}

uint64_type &uint64_type::operator=(BTProcessable &&other)
{
    val = reinterpret_cast<uint64_type &&>(other).val;

    return *this;
}

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

    uint64_type initobj{0};

    B_Tree tr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 0);
    Node nd1(&initobj, tr.maxNumOfObjectsInNode, 45);
    tr.initNode(nd1);

    nd1._numOfCurrentStoredObjects = nd1.nodesValPtrs.size();

    for (size_t i = 0; i < nd1.nodesValPtrs.size(); i++)
        nd1.nodesValPtrs[i] = new uint64_type(i * 9487 + 234);

    for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
        nd1.childrenNodesIds[i] = (i * 9875 + 34);

    tr._cache.saveNode(nd1);

    Node nd2(&initobj, tr.maxNumOfObjectsInNode, 0);
    tr._cache.getNode(2, nd2);

    ASSERT_EQ(nd1.childrenNodesIds.size(), nd2.childrenNodesIds.size());
    ASSERT_EQ(nd1.nodesValPtrs.size(), nd2.nodesValPtrs.size());

    EXPECT_EQ(nd1.isLeaf, nd2.isLeaf);
    EXPECT_EQ(nd1._id, nd2._id);
    EXPECT_EQ(nd1._numOfCurrentStoredObjects, nd2._numOfCurrentStoredObjects);
    EXPECT_EQ(nd1._levelMarker, nd2._levelMarker);

    for (size_t i = 0; i < nd1.nodesValPtrs.size(); i++)
        EXPECT_EQ(*nd1.nodesValPtrs[i], *nd2.nodesValPtrs[i]);

    for (size_t i = 0; i < nd1.childrenNodesIds.size(); i++)
        EXPECT_EQ(nd1.childrenNodesIds[i], nd2.childrenNodesIds[i]);
}

TEST(BtreeTests, TestSomeValues)
{
    size_t pageSizeInBytes = 64;

    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 0);

    std::vector<uint64_t> vals = {9, 10, 3, 134, 104, 17};

    for (auto val : vals)
    {
        uint64_type temp{val};
        btr.insertKey(temp);
    }

    for (auto val : vals)
    {
        uint64_type temp{val};
        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

        uint16_t index = 0;
        btr.search(btr._root, temp, findNode, index);

        ASSERT_NE(index, 0);
        ASSERT_EQ(*findNode.nodesValPtrs[index - 1], temp);
    }
}

TEST(BtreeTests, TestRandomValues)
{
    srand(std::time(NULL));
    std::vector<uint64_t> vals;

    createRandomShuffle(1, 100'000, 40, vals);

    size_t pageSizeInBytes = 128;
    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 4);

    for (auto val : vals)
    {
        uint64_type temp{val};
        btr.insertKey(temp);
    }

    for (auto val : vals)
    {
        uint64_type temp{val};
        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

        uint16_t index = 0;
        btr.search(btr._root, temp, findNode, index);

        ASSERT_NE(index, 0);
        ASSERT_EQ(*findNode.nodesValPtrs[index - 1], temp);
    }
}

TEST(BtreeTests, TestMoreValues)
{
    srand(std::time(NULL));
    std::vector<uint64_t> sh;

    createRandomShuffle(1, 100'000, 2000, sh);

    size_t pageSizeInBytes = 4096;
    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 4);

    for (auto i : sh)
    {
        uint64_type temp{i};
        btr.insertKey(temp);
    }

    for (auto i : sh)
    {
        uint64_type temp{i};
        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

        uint16_t index = 0;
        btr.search(btr._root, temp, findNode, index);

        ASSERT_NE(index, 0);
        ASSERT_EQ(*findNode.nodesValPtrs[index - 1], temp);
    }
}

TEST(BtreeTests, TestDeletionOfSomeValues)
{
    size_t pageSizeInBytes = 64;
    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 0);

    std::vector<uint64_type> vals = {9, 10, 3, 134, 104, 17};

    for (auto val : vals)
        btr.insertKey(val);

    for (auto val : vals)
    {
        btr.deleteKey(btr._root, val);

        uint64_type searchVal = val;
        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

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
    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, true, std::string("nodes"), pageSizeInBytes, 4);

    for (auto val : vals)
    {
        uint64_type temp{val};
        btr.insertKey(temp);
    }

    for (auto val : vals)
    {
        uint64_type temp{val};
        btr.deleteKey(btr._root, temp);

        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

        uint16_t index = 0;
        btr.search(btr._root, temp, findNode, index);

        ASSERT_EQ(index, 0);
    }
}

TEST(BtreeTests, TestDeletionOfMoreValues)
{
    srand(std::time(NULL));
    std::vector<uint64_t> vals;

    createRandomShuffle(1, 1'000'000, 15'000, vals);

    size_t pageSizeInBytes = 4096;
    uint64_type initobj{0};

    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 5);

    for (auto val : vals)
    {
        uint64_type temp{val};
        btr.insertKey(temp);
    }

    for (auto val : vals)
    {
        uint64_type temp{val};

        btr.deleteKey(btr._root, temp);

        Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

        uint16_t index = 0;
        btr.search(btr._root, temp, findNode, index);

        ASSERT_EQ(index, 0);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}