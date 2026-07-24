#include "lcg2.h"
#include "Parsing.h"
#include "Assembly.h"
#include "2D7G.h"

int main() {
	Manager manager;
	
	manager.initialize(10000);
	manager.iteration();
	return 0;
}

