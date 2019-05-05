#include "Agent.h" 

using namespace std;

Agent::Agent(void)
	:dim(0), fit(0)
{

}

Agent::~Agent(void)
{
}

void Agent::initialize(int d)
{
	dim = d;
	fit = 0;
	if (dim > 0)
	{
		pos.reserve(d);
		pos = random(down, up, dim);
	}
}
