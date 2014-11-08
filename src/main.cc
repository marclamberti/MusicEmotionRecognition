#include <aubio/aubio.h>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace {
	std::vector<std::vector<float>> features;

	const int size = 512;
}


/**
 * Read the sound file and display a block of hop_size samples
 * @param  file the soung file
 * @return      true if there was no error, otherwise false
 */
bool	readWav(const std::string &file) {
	uint_t 	samplerate = 0;
  	uint_t 	hop_size = 256;
  	uint_t 	n_frames = 0;
  	uint_t	read = 0;

  	aubio_source_t *aubio_source = new_aubio_source(const_cast<char *>(file.c_str()), samplerate, hop_size);
  	if (!aubio_source) {
  		std::cerr << "new_aubio_source failed" << std::endl;
  		return false;
  	}
  	
  	fvec_t *vec = new_fvec(hop_size);
  	samplerate =  aubio_source_get_samplerate(aubio_source);

  	do {
  		aubio_source_do(aubio_source, vec, &read);
  		fvec_print(vec);
  		n_frames += read;
  	} while (read == hop_size);

  	std::cout << "read " << n_frames << " frames at " << samplerate << "Hz (" << n_frames / hop_size << " blocks) from " << file << std::endl;
	del_fvec(vec);
 	del_aubio_source(aubio_source);
  	return true;
}


std::string	GetFileExtension(std::string file) {
	std::string::size_type idx = file.rfind('.');
	if (idx != std::string::npos)
		return file.substr(idx + 1);;
	return std::string("");
}

int kNumberOfJazzSong = 85;
int kNumberOfClassicalSong = 65;

enum Format {
	MIDI,
	WAV,
	NUMBER_OF_FORMAT,
};

std::vector<std::string> kFormatToString = {
	".mid",
	".wav",
};

enum Genders {
	CLASSICAL,
	JAZZ,
	NUMBER_OF_GENDERS
};

enum Modes {
	MINOR,
	MAJOR,
	NUMBER_OF_MODES
};

enum Keys {
	C,
	DB,
	D,
	EB,
	E,
	F,
	GB,
	G,
	AB,
	A,
	BB,
	B,
	NUMBER_OF_KEYS
};

std::map<std::string, enum Keys> kStringToKey = {
	{"C", Keys::C},
	{"DB", Keys::DB},
	{"D", Keys::D},
	{"EB", Keys::EB},
	{"E", Keys::E},
	{"F", Keys::F},
	{"GB", Keys::GB},
	{"G", Keys::G},
	{"AB", Keys::AB},
	{"A", Keys::A},
	{"BB", Keys::BB},
	{"B", Keys::B}
};

int kModeAndKeyIndex = 10;
int kTempoIndex = 6;

enum Features {
	BPM,
	AVERAGE_ENERGY,
	ENERGY_STANDARD_DEVIATION,
	AVERAGE_FUNDAMENTAL_FREQUENCY,
	FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION,
	NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY,
	AVERAGE_CENTROID,
	CENTROID_STANDARD_DEVIATION,
	MODE,
	KEY,
	GENDER,
	NUMBER_OF_FEATURES
};

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> Split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string 	ExecCommand(std::string const &cmd) {
	int		status;
	pid_t	pid;
	std::string output;
	
	FILE *lsofFile_p = popen(cmd.c_str(), "r");
  	if (lsofFile_p)
  	{  
  		char buffer[1024];
  		memset(buffer, 0, 1024);
  		// For now we skip the first line.
  		char *line_p = fgets(buffer, sizeof(buffer), lsofFile_p);
  		line_p = fgets(buffer, sizeof(buffer), lsofFile_p);
  		pclose(lsofFile_p);
  		return line_p;
  	}

	return output;
}

std::string	IdToSongName(int song_id, enum Format f) {
	std::string song_number;
	std::string song_genre;

	if (song_id <= kNumberOfClassicalSong) {
		song_genre = "clas";
		if (song_id >= 10)
			song_number = std::to_string(song_id);
		else
			song_number = "0" + std::to_string(song_id);
	} else {
		song_genre = "jazz";
		if (song_id - kNumberOfClassicalSong >= 10)
			song_number = std::to_string(song_id - kNumberOfClassicalSong);
		else
			song_number = "0" + std::to_string(song_id - kNumberOfClassicalSong);
	}
	return song_genre + "_" + song_number + kFormatToString[f];
}

std::string GetInfoFromMIDI(std::string const &song_path) {
	char buff[512];
	memset(buff, 0, 512);
	getcwd(buff, 512);

	std::string cmd = "metamidi -l ";
	std::string cwd = buff;
	return ExecCommand(cmd + cwd + "/" + song_path);
}

enum Genders ExtractGender(int song_id) {
	return song_id <= kNumberOfClassicalSong ? Genders::CLASSICAL : Genders::JAZZ;
}

float ExtractBpm(std::vector<std::string> const &midi_info) {
	return std::atof(midi_info[kTempoIndex].c_str());
}

enum Modes ExtractMode(std::vector<std::string> const &midi_info) {
	if (midi_info[kModeAndKeyIndex].find('m') != std::string::npos)
		return Modes::MINOR;
	return Modes::MAJOR;
}

enum Keys ExtractKey(std::vector<std::string> const &midi_info) {
	std::string info_line = midi_info[kModeAndKeyIndex];
	std::transform(info_line.begin(), info_line.end(), info_line.begin(), ::toupper);
	std::vector<std::string> key_and_mode = Split(info_line, 'M');

	for (auto map_element : kStringToKey) {
		if (key_and_mode[0] == map_element.first)
			return map_element.second;
	}
	return Keys::NUMBER_OF_KEYS;
}

void	FillFeatures(std::vector<float> &tuple, int song_id, std::string const &song_path) {
	std::cout << "Name is : " << song_path << std::endl;
	tuple.resize(NUMBER_OF_FEATURES + 2);
	std::vector<std::string> info_in_midi = Split(GetInfoFromMIDI(song_path), ';');
	tuple[BPM] = ExtractBpm(info_in_midi);
	tuple[MODE] = static_cast<float>(ExtractMode(info_in_midi));
	tuple[KEY] = static_cast<float>(ExtractKey(info_in_midi));
	tuple[GENDER] = static_cast<float>(ExtractGender(song_id));

	std::cout << "The beat is : " << tuple[BPM] << std::endl;
	std::cout << "The mode is : " << tuple[MODE] << std::endl;
	std::cout << "The key is : " << tuple[KEY] << std::endl;
	std::cout << "The gender is : " << tuple[GENDER] << std::endl;
}

int kValenceIndex = 11;
int kArousalIndex = 12;
int kLabelIndex = std::min(kArousalIndex, kValenceIndex);

enum LabelTypes {
	AROUSAL,
	VALENCE,
	MULTILABEL,
	NUMBER_OF_LABEL_TYPES,
};

std::string FormatTuple(std::vector<float> const &tuple, enum LabelTypes lt) {
	std::string formatted_tuple;

	if (lt == LabelTypes::AROUSAL || lt == LabelTypes::MULTILABEL)
		formatted_tuple += std::to_string(tuple[kArousalIndex]);
	if (lt == LabelTypes::VALENCE || lt == LabelTypes::MULTILABEL)
		formatted_tuple += std::to_string(tuple[kArousalIndex]);

	int i = 1;
	for (auto const &attribute : tuple) {
		if (i <= kLabelIndex)
			formatted_tuple += " " + std::to_string(i++) + ':' + std::to_string(attribute);
	}
	return formatted_tuple;
}

void FormatDatasetAndWriteInFile(std::vector<std::vector<float>> const &data_set, std::string const &filename,
								 enum LabelTypes lt) {
	std::ofstream output_file_stream(filename.c_str(), std::fstream::out);

	if (!output_file_stream.is_open()) {
		std::cerr << "Cannot open the file " << filename << std::endl;
		return;
	}

	for (auto const &tuple : data_set) {
		std::string line = FormatTuple(tuple, lt);
		std::cout << "The line : " << line << std::endl;
		line += '\n';
		output_file_stream << line;
	}
	output_file_stream.close();
}

int main(int ac, char **av){
	std::vector<std::vector<float>> data_set(kNumberOfClassicalSong + kNumberOfJazzSong); 
	
	if (ac < 2) {
		std::cout << "Usage:" << av[0] << " <training_directory>" << std::endl;
		return 1;
	}
	int song_id = 1;
	for (auto &tuple : data_set) {
		std::string song_path = av[1] + std::string("/") + IdToSongName(song_id, Format::MIDI);
		FillFeatures(tuple, song_id++, song_path);
	}
	FormatDatasetAndWriteInFile(data_set, "test.txt", LabelTypes::AROUSAL);
    return 0;
}