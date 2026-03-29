#include "btprocessable.h"

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