#include "lcg2.h"
#include "Parsing.h"
#include "Assembly.h"
#include "2D7G.h"

int main() {
	Manager manager;

	manager.initialize();
	manager.iteration();
	return 0;
}

