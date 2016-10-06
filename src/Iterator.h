#ifndef __ITERATOR_H__
#define __ITERATOR_H__

#include <vector>
#include <stdlib.h>

template <typename Item>
class Iterator
{
public:
	Iterator()
	{

	}

	virtual ~Iterator()
	{

	}

public:
	virtual Item* GetFirstItem() = 0;
	virtual Item* GetNextItem() = 0;
	virtual Item* GetCurrentItem() = 0;
};

template<typename Item>  
class ConcreteAggregate;  

template<typename Item>  
class ConcreteIterator : public Iterator <Item>  
{  
public:  
	ConcreteIterator(ConcreteAggregate<Item>* a)
		: aggr(a), cur(0)
	{

	}  

	virtual Item* GetFirstItem()  
	{  
		cur = 0;
		return &(*aggr)[0];
	}  

	virtual Item* GetNextItem()  
	{  
		if(++cur < aggr->GetLen())
		{
			return &(*aggr)[cur];
		}
		else
		{
			return NULL;
		}
	}

	virtual Item* GetCurrentItem()
	{
		if(cur < aggr->GetLen())
			return &(*aggr)[cur];
		else
			return NULL;
	}
private:
	ConcreteAggregate<Item> * aggr;  
	int cur;  
};  

template<typename Item>  
class Aggregate  
{  
public:  
	virtual Iterator<Item>* CreateIterator() = 0;  
	virtual ~Aggregate()
	{

	}  
};  

template<typename Item>  
class ConcreteAggregate : public Aggregate<Item>  
{  
public:  
	ConcreteAggregate(const std::vector<Item>& data)  
	{  
		m_data = data;
	}  

	virtual Iterator<Item>* CreateIterator()  
	{  
		return new ConcreteIterator<Item>(this);  
	}  

	Item& operator[](const int index)
	{  
		return m_data[index];  
	}  

	int GetLen() const
	{  
		return m_data.size();  
	}  

private:
	std::vector<Item> m_data;
};  

#endif // __ITERATOR_H__

