#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <vector>

std::vector<std::string> getFileList(std::string filePath) {
    std::vector<std::string> files;
    std::ifstream inputFile(filePath.c_str());
    std::string line;
    while (std::getline(inputFile, line)) {
        files.push_back(line);
    }
    inputFile.close();
    return files;
}

std::map<std::string, int> countKeys(std::vector<std::string> files) {
    std::map<std::string, int> keyCounts;
    for (int i = 0; i < files.size(); i++) {
        std::ifstream file(files[i].c_str());
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key;
            int value;
            if (!(iss >> key >> value)) {
                std::cout << "Error reading line: " << line << std::endl;
                break;
            } 
            keyCounts[key] += value;
        }
        file.close();
    }
    return keyCounts;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " [-f file_list.txt | file1.txt file2.txt ...]\n";
        return 1;
    }

    std::vector<std::string> files;
    if (std::string(argv[1]) == "-f") {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " -f file_list.txt\n";
            return 1;
        }
        files = getFileList(argv[2]);
    } else {
        files.assign(argv + 1, argv + argc);
    }

    std::map<std::string, int> keyCounts = countKeys(files);

    for (const auto& pair : keyCounts) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }

    return 0;
}
