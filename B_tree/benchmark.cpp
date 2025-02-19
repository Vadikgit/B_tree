#include <benchmark/benchmark.h>
#include "btree.h"
#include <vector>
#include <thread>

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

void InsertionOfValues(B_Tree<uint64_t> &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
{
  srand(randomSeed);
  std::vector<uint64_t> vals;

  std::vector<uint64_t> partitionIndexesShuffle;
  createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

  for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // insertion
  {
    createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

    for (auto val : vals)
      btr.insertKey(val);
  }
}

void SearchOfValues(B_Tree<uint64_t> &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
{
  srand(randomSeed);
  std::vector<uint64_t> vals;

  std::vector<uint64_t> partitionIndexesShuffle;
  createRandomShuffle(0, numOfValuesToProcess / numOfValuesInOneSubShuffle - 1, numOfValuesToProcess / numOfValuesInOneSubShuffle, partitionIndexesShuffle);

  for (size_t i = 0; i < partitionIndexesShuffle.size(); i++) // search
  {
    createRandomShuffle(partitionIndexesShuffle[i] * numOfValuesInOneSubShuffle, (partitionIndexesShuffle[i] + 1) * numOfValuesInOneSubShuffle - 1, numOfValuesInOneSubShuffle, vals);

    for (auto val : vals)
    {
      uint64_t searchVal = val;
      Node<uint64_t> findNode(btr.maxNumOfObjectsInNode, 0);

      uint16_t index = 0;
      btr.search(btr._root, searchVal, findNode, index);
    }
  }
}

void DeletionOfValues(B_Tree<uint64_t> &btr, size_t numOfValuesToProcess, size_t numOfValuesInOneSubShuffle, unsigned int randomSeed)
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
      btr.deleteKey(btr._root, val);
    }
  }
}

static void BM_BTreeInsertion(benchmark::State &state)
{
  size_t numOfValuesToProcess = state.range(0);
  size_t numOfValuesInOneSubShuffle = state.range(1);

  unsigned int randSeed = 0;
  for (auto _ : state)
  {
    size_t pageSizeInBytes = 4096;
    B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
    InsertionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, randSeed++);
  }
}

static void BM_BTreeSearch(benchmark::State &state)
{
  size_t numOfValuesToProcess = state.range(0);
  size_t numOfValuesInOneSubShuffle = state.range(1);

  size_t pageSizeInBytes = 4096;
  B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
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

  unsigned int randSeed = 0;
  for (auto _ : state)
  {
    state.PauseTiming();
    size_t pageSizeInBytes = 4096;
    B_Tree<uint64_t> btr(false, false, std::string("nodes"), pageSizeInBytes, 1 + 250 + 250 * 250);
    InsertionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, 0);
    state.ResumeTiming();

    DeletionOfValues(btr, numOfValuesToProcess, numOfValuesInOneSubShuffle, randSeed++);
  }
}

BENCHMARK(BM_BTreeInsertion)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});
BENCHMARK(BM_BTreeSearch)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});
BENCHMARK(BM_BTreeDeletion)->ArgsProduct({benchmark::CreateRange(10'000, 1'000'000, 10), {2000}});

BENCHMARK_MAIN();