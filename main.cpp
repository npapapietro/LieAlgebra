#include"libs/RootStructure.h"

/*
* This is a test file to ensure all headers are correctly set on your system.
*/
using namespace std;
using namespace Representation;

int main(int argc, char const ** argv)
{

	if(argc > 1){
		std::string loglevel = argv[1];
		Diagnositics::Log::ReportingLevel() = Diagnositics::Log::FromString(loglevel);
	}
	else{
		Diagnositics::Log::ReportingLevel() = Diagnositics::logINFO;

	}

    LieBase<GroupType::G> f(2);

    // for (const auto& i: f.get_positiver())
    // {
    //     cout<<i.ortho.transpose()<<endl;
    // }
    cout<<f.get_Cartan()<<endl;
    
    return 0;
}