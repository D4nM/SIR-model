set -e
g++ SIR_model_v2.cpp -c -Wextra -Wall -fsanitize=address -fno-omit-frame-pointer
g++ SIR_model_v2.o -lsfml-graphics -lsfml-window -lsfml-system -fsanitize=address -o SIR_model_v2_executable
./SIR_model_v2_executable
