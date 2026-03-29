#pragma once

#include <vector>
#include <functional>
#include <memory>
#include <cstdint>

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

static std::function<bool(const std::unique_ptr<BTProcessable> &, const std::unique_ptr<BTProcessable> &)> BTProcessablePtrComp = [](const std::unique_ptr<BTProcessable> &el1, const std::unique_ptr<BTProcessable> &el2)
{
    return (*el1 < *el2);
};
