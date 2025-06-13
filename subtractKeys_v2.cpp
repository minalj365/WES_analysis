#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <iostream>
#include <vector>

std::map<std::string, int> loadCounts(std::string filename) {
    std::map<std::string, int> keyCounts;
    std::ifstream file(filename.c_str());
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        std::string value_str;
        if (std::getline(iss, key, ':')) {
            iss.ignore(); // Ignore the space after the colon
            if (!(iss >> value_str)) {
                std::cout << "Error reading line: " << line << std::endl;
                break;
            }
            int value;
            try
            {
                value = std::stoi(value_str); // Convert string to int
            }
            catch(const std::invalid_argument& e)
            {
                 std::cerr << "Invalid argument: " << value_str << " in line: " << line << std::endl;
                 continue;
            } catch (const std::out_of_range& e) {
                 std::cerr << "Out of range: " << value_str << " in line: " << line << std::endl;
                 continue;
            }
              
            keyCounts[key] = value;
        }
    }
    file.close();
    return keyCounts;
}

void subtractKeys(std::string filename, std::map<std::string, int>& keyCounts) {
    std::ifstream file(filename.c_str());
    std::string key;
    while (std::getline(file, key)) {
        keyCounts.erase(key);
    }
    file.close();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <counts_file.txt> <keys_file1.txt> <keys_file2.txt> ...\n";
        return 1;
    }

    std::map<std::string, int> keyCounts = loadCounts(argv[1]);

    // Subtract keys from all keys files
    for (int i = 2; i < argc; i++) {
        subtractKeys(argv[i], keyCounts);
    }

    for (const auto& pair : keyCounts) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }

    return 0;
}
