#define _CRT_SECURE_NO_WARNINGS
#include "btree.h"

struct point : public BTProcessable
{
	uint64_t x{0};
	uint64_t y{0};

	virtual void fromBytes(const std::vector<uint8_t> &ibuf, size_t firstBytePos) override;
	virtual void toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos) override;
	virtual bool less(const BTProcessable &w2) const override;
	virtual size_t sizeInBytes() override;
	virtual point *createNew() override;

	point() : BTProcessable() {};
	point(uint64_t xx, uint64_t yy) : x{xx}, y{yy}, BTProcessable() {};
	point(const point &other) : x{other.x}, y{other.y}, BTProcessable(other) {};
	point(point &&other) : x{other.x}, y{other.y}, BTProcessable(other) {};

	virtual point &operator=(const point &other);
	virtual point &operator=(point &&other);

	virtual point &operator=(const BTProcessable &other) override;
	virtual point &operator=(BTProcessable &&other) override;
	virtual ~point() {};
};

void point::fromBytes(const std::vector<uint8_t> &ibuf, size_t firstBytePos)
{
	x = *(reinterpret_cast<const uint64_t *>(&ibuf[firstBytePos]));
	firstBytePos += sizeof(uint64_t);
	y = *(reinterpret_cast<const uint64_t *>(&ibuf[firstBytePos]));
}

void point::toBytes(std::vector<uint8_t> &obuf, size_t firstBytePos)
{
	*(reinterpret_cast<uint64_t *>(&obuf[firstBytePos])) = x;
	firstBytePos += sizeof(uint64_t);
	*(reinterpret_cast<uint64_t *>(&obuf[firstBytePos])) = y;
}

bool point::less(const BTProcessable &p2) const
{
	if (x != (dynamic_cast<const point &>(p2)).x)
		return x < (dynamic_cast<const point &>(p2)).x;

	return y < (dynamic_cast<const point &>(p2)).y;
}

size_t point::sizeInBytes()
{
	return sizeof(x) + sizeof(y);
}

point *point::createNew()
{
	return new point();
}

point &point::operator=(const point &other)
{
	if (this != &other)
	{
		x = other.x;
		y = other.y;
	}

	return *this;
}

point &point::operator=(point &&other)
{
	x = other.x;
	y = other.y;

	return *this;
}

point &point::operator=(const BTProcessable &other)
{
	if (this != &other)
	{
		x = reinterpret_cast<const point &>(other).x;
		y = reinterpret_cast<const point &>(other).y;
	}

	return *this;
}

point &point::operator=(BTProcessable &&other)
{
	x = reinterpret_cast<point &&>(other).x;
	y = reinterpret_cast<point &&>(other).y;

	return *this;
}

int main()
{
	point tp1(1, 2);

	size_t pageSizeInBytes = 128;
	B_Tree tree{&tp1, tp1.sizeInBytes(), false, false, "nodes", pageSizeInBytes, 4};

	tree.insertKey(tp1);

	auto searchVal = tp1;
	Node findNode(&tp1, tree.maxNumOfObjectsInNode, 0);

	uint16_t index = 0;
	tree.search(tree._root, searchVal, findNode, index);

	std::cout << index << " " << findNode.nodesValPtrs[index - 1] << " " << &tp1 << std::endl;
	std::cout << " " << (*dynamic_cast<point *>(findNode.nodesValPtrs[index - 1])).x << " " << (*dynamic_cast<point *>(findNode.nodesValPtrs[index - 1])).y << std::endl;

	return 0;
}
