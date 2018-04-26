#ifndef TMA_TYPES_H_
#define TMA_TYPES_H_

namespace tma
{

typedef unsigned int uint;
typedef double       real;



template<typename key, typename value> 
class map : public std::map<key, value>	// own map class that extends std::map
{
public:
	bool contains(value v) const
	{
		for(auto it = this->begin(); it != this->end(); ++it) {
			if(it->second == v) 
				return true;
		}
		return false;
	}
	
	friend std::ostream& operator <<(std::ostream& os, const map& m)
	{
		os << "[";
		for(auto it = m.begin(); it != m.end(); ++it)
			os << it->first << ":" << it->second << ( it != --m.end() ? "; " : "]");
		return os;
	}
};

template<typename T>
class vector : public std::vector<T>
{
public:
	bool contains(T v) const
	{
		for(uint i = 0; i < this->size(); i++) {
			if(this->at(i) == v)
				return true;
		}
		return false;
	}
	
	uint maxIndex() const
	{
		uint m = 0;
		for(uint i = 1; i < this->size(); i++) {
			if(this->at(i) > this->at(m))
				m = i;
		}
		return m;
	}
	
	uint find(T v) const
	{
		for(uint i = 0; i < this->size(); i++) {
			if(this->at(i) == v)
				return (int)i;
		}
		return this->size();
	}
	
	using std::vector<T>::vector;
	
	friend std::ostream& operator <<(std::ostream& os, const vector& v)
	{
		os << "[";
		for(auto it = v.begin(); it != v.end(); ++it)
			os << (*it) <<  ( it != --v.end() ? "\t" : "]");
		return os;
	}
};

} /* namespace tma */

#endif /* TMA_TYPES_H _ */
