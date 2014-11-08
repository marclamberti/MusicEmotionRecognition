#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

namespace {
	int kNumberOfJazzSong = 85;
	int kNumberOfClassicalSong = 65;
	int kTotalData = kNumberOfJazzSong + kNumberOfClassicalSong;

	enum Intension {
		READ,
		WRITE,
	};

	std::map<enum Intension, std::ios_base::openmode> kIntensionToIFstream { 
		{READ, std::ifstream::in},
		{WRITE, std::ifstream::out},
	};
}

bool Open(char *filename, enum Intension intension, std::ifstream &fs) {
	fs.open(filename, kIntensionToIFstream[intension]);
	return fs.is_open();
}

std::pair<float, float> GetValues(std::string const &line) {
	std::pair<float, float> values;
	std::istringstream iss(line);
	std::string tmp;

	iss >> tmp;
	values.first = std::atof(tmp.c_str());
	iss >> tmp;
	values.second = std::atof(tmp.c_str());
	return values;
}

int main(int ac, char **av) {
	if (ac < 5) {
		std::cerr << 	"usage : ./merge_results <filename Jordy> <filename"
						" Marc> <filename Rossi> <output file>" << std::endl;
	}
	std::vector<std::ifstream> input_file_streams(3);
	std::ofstream output_file_stream(av[ac - 1], std::fstream::out);

	if (!output_file_stream.is_open()) {
		std::cerr << "Cannot open the file " << av[ac - 1] << std::endl;
		return 1;
	}

	for (int i = 1; i < ac - 1; ++i) {
		if (!Open(av[i], Intension::READ, input_file_streams[i - 1])) {
			std::cerr << "Cannot open the file " << av[i] << std::endl;
			return 1;
		}
	}

	std::vector<std::pair<float, float>> final_output(kTotalData);

	for (auto &tuple : final_output) {
		for (auto &file : input_file_streams) {
			std::string line;
			file >> line;
			std::cout << "The line is : " << line << std::endl;
			std::pair<float, float> values = GetValues(line);
			tuple.first += values.first;
			tuple.second += values.second;
		}
		tuple.first /= input_file_streams.size();
		tuple.second /= input_file_streams.size();
	}

	for (auto &values : final_output) {
		std::string line = std::to_string(values.first) + " " + std::to_string(values.second); 
		output_file_stream << line + '\n';
		std::cout << "line : " <<  line << std::endl;
	}

	for (auto &file : input_file_streams) {
		file.close();
	}

	output_file_stream.close();
}