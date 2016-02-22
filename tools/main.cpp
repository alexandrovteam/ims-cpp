#include <iostream>
#include <map>
#include <string>

typedef std::map<std::string, int(*)(int, char**)> SubCmdMap;
#define _(x) extern int x##_main(int, char**)
#define __(x) {#x, &x##_main}

 _(convert);  _(isocalc);  _(detect);

const SubCmdMap subcommands{
__(convert), __(isocalc), __(detect)
};

#undef _
#undef __

int main(int argc, char** argv) {
	if (argc == 1) {
		std::cout << "Usage: ims <subcommand> [options...]\n"
							<< "The following subcommands are available:\n"
							<< "    convert - converts .imzML to .imzb format\n"
							<< "    isocalc - computes isotope patterns for a list of formulas\n"
							<< "    detect  - scores metabolite-spectrum matches for a dataset\n"
							<< "\n"
							<< "To get help for a subcommand, run it without any options."
							<< std::endl;
		return 0;
	}

	auto it = subcommands.find(argv[1]);
	if (it == subcommands.end()) {
		std::cout << "Unknown subcommand '" << argv[1] << "'" << std::endl;
		return -1;
	}

	return it->second(argc - 1, argv + 1);
}
