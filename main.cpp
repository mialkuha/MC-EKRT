// A Monte Carlo EKRT event generator, see https://doi.org/10.48550/arXiv.2406.17592
// <Citation>.
// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.
//
// This program uses OpenMP, which is licensed under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "mcaa.hpp"

int main(int argc, char *argv[])
{
    std::cout << "MC-EKRT: Copyright (c) 2025 Mikko Kuha (University of Jyväskylä)." << std::endl;
    std::cout << "This program is published under GPLv3 license, and comes with" << std::endl; 
    std::cout << "ABSOLUTELY NO WARRANTY; This is free software, and you are welcome" << std::endl;
    std::cout << "to redistribute it under certain conditions; For more information," << std::endl;
    std::cout << "see <http://www.gnu.org/licenses/>" << std::endl;

    if (argc < 2)
    {
        mcaa sim("params_template");
        sim.run();
    }
    else
    {
        const std::string filename{argv[1]};
        mcaa sim(filename);
        sim.run();
    }

    return 0;
}
