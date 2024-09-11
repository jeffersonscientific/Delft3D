#pragma once
#include <list>

struct keyValueLL {
	char * key;
	char * val;
	keyValueLL * nextkv;
};
struct keyValue {
	char * key;
	char * val;
};

typedef std::list<keyValue>	keyValueList;
