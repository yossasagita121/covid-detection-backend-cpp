#include <iostream>
#include <iomanip>
#include "json.hpp"
#include "httplib.h"
#include "functions.h"
#include <iostream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

using namespace httplib;
using json = nlohmann::json;
using namespace std;

string exec(string cmd)
{
    unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(cmd.c_str(), "r"), _pclose);
    if (!pipe)
    {
        throw runtime_error("popen() failed!");
    }
    string result;
    char buffer[128];
    while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr)
    {
        result += buffer;
    }
    return result;
}

int main()
{
    // Membuat server
    Server svr;

    // Menangani route GET di root "/"
    svr.Get("/", [](const Request &req, Response &res)
            { res.set_content("Hello, World!", "text/plain"); });

    svr.Post("/run", [](const Request &req, Response &res)
             { 
        string upload_dir = "dataseries/";
        auto file = req.files.find("file");

        if (file != req.files.end()) {
            const httplib::MultipartFormData& file_data = file->second;
            ofstream output_file(upload_dir + file_data.filename, ios::binary);
            cout << file_data.filename << endl;
            
            if (output_file) {
                output_file.write(file_data.content.c_str(), file_data.content.size());
                output_file.close();

                string command = "main.exe dataseries/" + file_data.filename;

                try
                {
                    string output = exec(command);
                    res.set_content(output, "text/plain");
                }
                catch (const exception &e)
                {
                    cerr << "Error: " << e.what() << endl;
                }
            } else {
                res.status = 500;
                res.set_content("Failed to save file", "text/plain");
            }
        } else {
            res.status = 400;
            res.set_content("No file uploaded", "text/plain");
        } });

    // Menangani route POST untuk API "/api/data"
    svr.Post("/data", [](const Request &req, Response &res)
             {
        string received_data = req.body;
        string response_data = received_data;

        // Parse string JSON menjadi objek JSON
        json j = json::parse(received_data);

        // Mendapatkan ukuran dimensi array
        int rows = j["data"].size();
        int cols = j["data"][0].size();

        // Menginisialisasi array 2D
        double data[DIMENSION][SAMPLE_MAX];

        // Mengisi array 2D dengan data dari JSON
        for (int i = 0; i < rows; ++i)
        {
            for (int k = 0; k < cols; ++k)
            {
                data[i][k] = j["data"][i][k];
            }
        }

        // Menampilkan 5 data pertama dari array 2D untuk memverifikasi
        int count = 0;
        for (int i = 0; i < rows && count < 5; ++i)
        {
            for (int k = 0; k < cols && count < 5; ++k)
            {
                cout << fixed << setprecision(10) << "[" << data[i][k] << "]" << endl;
                ++count; // Tambahkan count setiap kali data ditampilkan
            }
        }

        string result = native_lib_main(data);

        res.set_content(result, "text/plain"); });

    svr.Post("/file", [](const Request &req, Response &res)
             {
        string received_data = req.body;
        cout << "File Content = " << received_data << endl; });

    // Tampilkan pesan bahwa server sedang berjalan
    cout << "Server running on port 8080" << endl;

    // Jalankan server di port 8080
    svr.listen("0.0.0.0", 8080);

    return 0;
}
