#include <benchmark/benchmark.h>
#include "../btree.h"
#include <vector>
#include <thread>

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

void InsertionOfValues(B_Tree &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
{
    srand(randomSeed);
    std::vector<uint64_t> vals;

    std::vector<uint64_t> partitionIndexesShuffle;
    createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

    for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // insertion
    {
        createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

        for (auto val : vals)
        {
            uint64_type temp{val};
            btr.insertKey(temp);
        }
    }
}

void SearchOfValues(B_Tree &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
{
    srand(randomSeed);
    std::vector<uint64_t> vals;

    uint64_type initobj{0};

    std::vector<uint64_t> partitionIndexesShuffle;
    createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

    for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search
    {
        createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

        for (auto val : vals)
        {
            uint64_type temp{val};
            Node findNode(&initobj, btr.maxNumOfObjectsInNode, 0);

            uint16_t index = 0;
            btr.search(btr._root, temp, findNode, index);
        }
    }
}

void DeletionOfValues(B_Tree &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
{
    srand(randomSeed);
    std::vector<uint64_t> vals;

    std::vector<uint64_t> partitionIndexesShuffle;
    createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

    for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // deletion
    {
        createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

        for (auto val : vals)
        {
            uint64_type temp{val};
            btr.deleteKey(btr._root, temp);
        }
    }
}

static void BM_BTreeInsertion(benchmark::State &state)
{
    size_t numOfValuesToProcess = state.range(0);
    size_t numOfValuesInOneSubShuffle = state.range(1);

    uint64_type initobj{0};

    unsigned int randSeed = 0;
    for (auto _ : state)
    {
        size_t pageSizeInBytes = 4096;
        B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
        InsertionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, randSeed++);
    }
}

static void BM_BTreeSearch(benchmark::State &state)
{
    size_t numOfValuesToProcess = state.range(0);
    size_t numOfValuesInOneSubShuffle = state.range(1);

    uint64_type initobj{0};

    size_t pageSizeInBytes = 4096;
    B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
    InsertionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, 0);

    unsigned int randSeed = 0;
    for (auto _ : state)
    {
        SearchOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, randSeed++);
    }
}

static void BM_BTreeDeletion(benchmark::State &state)
{
    size_t numOfValuesToProcess = state.range(0);
    size_t numOfValuesInOneSubShuffle = state.range(1);

    uint64_type initobj{0};

    unsigned int randSeed = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pageSizeInBytes = 4096;
        B_Tree btr(&initobj, initobj.sizeInBytes(), false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
        InsertionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, 0);
        state.ResumeTiming();

        DeletionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, randSeed++);
    }
}

BENCHMARK(BM_BTreeInsertion)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});
BENCHMARK(BM_BTreeSearch)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});
BENCHMARK(BM_BTreeDeletion)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});

BENCHMARK_MAIN();