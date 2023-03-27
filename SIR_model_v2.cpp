#include <iostream>
#include <vector>
#include <random> //for the random number generation
#include <cmath> //for math functions
#include <fstream>
#include <string>
#include <ctime>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <SFML/System.hpp> //for the clock in particular
#include <stdexcept>

constexpr double pi = 3.1415926536;

std::random_device rd; //seed
std::mt19937 gen{rd()};

enum class State
{
    Susceptible,
    Infected,
    Recovered
};

struct Position
{
    int R; //row
    int C; //column
};

bool operator==(Position const& leftPos, Position const& rightPos)
{
    return leftPos.R == rightPos.R && leftPos.C == rightPos.C;
}

template <typename T>
T gaussValue(float const mean, float const stddev)
{
    std::normal_distribution<T> gaussDist{mean, stddev};
    T result = gaussDist(gen);
    if(result < 0) { result = 0; }
    return result;
}

auto uniformIntValue(int const min, int const max)
{
    std::uniform_int_distribution<int> uniformDist{min, max};
    return uniformDist(gen);
}

auto uniformRealValue(float const min, float const max)
{
    std::uniform_real_distribution<float> uniformRealDist{min, max};
    return uniformRealDist(gen);
}

class Person
{
    private:
    State m_state = State::Susceptible;
    float m_resistance = 0.; //infection resistance -- ought to be a value within [0,1]
    float m_recoveryProbability = 0.; //ought to be a value within [0,1]
    bool m_inRange = false;
    bool m_hasMoved = false;

    public: //I'm not sure if delegating constructors could be used
    //Person(); // produces "class "Person" has no suitable copy constructor" error
    Person(); //defined this way I don't get the 'class "Person" has no suitable copy constructor' error
    
    //constructor with input fixed recovery probability
    Person(float const recov_probability) : Person(recov_probability, 0., 0.)
        {
        }

    //constructor with input resistance-generation parameters and fixed recovery probability
    Person(float const recov_probability,
           float const Res_mean,
           float const Res_variab) :
        m_resistance{gaussValue<float>(Res_mean, Res_variab)},
        m_recoveryProbability{recov_probability}
        {
            if(recov_probability < 0 || recov_probability >= 1)
            {
                throw std::runtime_error{"The recovery probability must be a probability, therefore a value between 0 and 1"};
            }
        }


    //constructor with input resistance-generation parameters and recovery probability generation parameters
    Person(float const recovProb_mean,
           float const recovProb_variab,
           float const Res_mean,
           float const Res_variab) :
        m_resistance{gaussValue<float>(Res_mean, Res_variab)},
        m_recoveryProbability{gaussValue<float>(recovProb_mean, recovProb_variab)}
        {}


    //State
    void set_state(State const& newState) { m_state = newState; }
    State const& get_state() const { return m_state; } //this method must be const (const after function name) in order for it to work in evolve()

    //Resistance
    void set_randomResistance(float const mean, float const stddev) //passed parameters must be within [0,1]
    {
        auto value = gaussValue<float>(mean, stddev);
        if(value < 0.) { value = 0.; }
        if(value > 1.) { value = 1.; }
        m_resistance = value;
    }
    auto get_resistance() const { return m_resistance; }

    //Recovery probability
    void set_recoveryProbability(float const probability)
    {
        m_recoveryProbability = probability;
    }
    void set_randomRecoveryProbability(float const mean, float const stddev)
    {
        auto value = gaussValue<float>(mean, stddev);
        if(value < 0.) { value = 0.; }
        if(value > 1.) { value = 1.; }
        m_recoveryProbability = value;
    }
    auto get_recoveryProbability() const { return m_recoveryProbability; }
    
    //inRange
    void set_inRange(bool const in) { m_inRange = in; }
    bool get_inRange() const { return m_inRange; }

    //hasMoved
    void set_hasMoved(bool const in) { m_hasMoved = in; }
    bool get_hasMoved() const { return m_hasMoved; }
    /*
    void display_attributes() const //for testing
    {
        std::string out;
        if(m_state == State::Susceptible) out = "Susceptible";
        if(m_state == State::Infected) out = "Infected";
        if(m_state == State::Recovered) out = "Recovered";
        std::cout << "State: " << out << '\n';
        std::cout << "Resist.: " << m_resistance << '\n';
        std::cout << "Recov. prob.: " << m_recoveryProbability << '\n';
        std::cout << "inRange: " << m_inRange << '\n';
        std::cout << "hasMoved: " << m_hasMoved << '\n';
    }
    */
};

class Board
{
    private:
    std::vector< std::vector<Person> > m_grid;    
    float m_infectionProbability = 0.;   
    int m_infectionRadius = 1;

    public:
    Board(int const length, int const height) : m_grid(height, std::vector<Person>(length, Person{}))
        {
            if(length <= 0 || height <= 0)
            {
                throw std::runtime_error{"Both grid dimensions must be a positive integer value"};
            }
        }

    Board(int const length, int const height, float const infect_prob, int const infect_radius) :
        m_grid(height, std::vector<Person>(length, Person{})),
        m_infectionProbability{infect_prob},
        m_infectionRadius{infect_radius}
        {
            if(length <= 0 || height <= 0)
            {
                throw std::runtime_error{"Both grid dimensions must be a positive integer value"};
            }
            if(infect_prob < 0 || infect_prob >= 1)
            {
                throw std::runtime_error{"The infection probability must be a probability, therefore a value between 0 and 1"};
            }
            if(infect_radius < 1)
            {
                throw std::runtime_error{"The infection radius must be a positive integer number"};
            }
        }
    
    Board(int const length, int const height, float const infect_prob, int const infect_radius, float const recov_prob) :
        m_grid(height, std::vector<Person>(length, Person{recov_prob})),
        m_infectionProbability{infect_prob},
        m_infectionRadius{infect_radius}
        {
            if(length <= 0 || height <= 0)
            {
                throw std::runtime_error{"Both grid dimensions must be a positive integer value"};
            }
            if(infect_prob < 0 || infect_prob >= 1)
            {
                throw std::runtime_error{"The infection probability must be a probability, therefore a value between 0 and 1"};
            }
            if(infect_radius < 1)
            {
                throw std::runtime_error{"The infection radius must be a positive integer number"};
            }
        }

    //Dimensions
    auto get_Length() const { return m_grid[0].size(); } // # of Columns, size of an 'inner' vector
    auto get_Heigth() const { return m_grid.size(); } // # of Rows, size of 'outer' vector

    //Grid
    auto const& get_Grid() const { return m_grid; }
    void set_Grid(std::vector< std::vector<Person> > const& newGrid) { m_grid = newGrid; }

    //Resistance of the population
    void setResistanceForAll(float const mean, float const stddev) //sets randomly the infection-resistance values of each Person
    {
        const int max_column = get_Length();
        const int max_row = get_Heigth();
        for(int row = 0; row < max_row; ++row)
        {
            for(int column = 0; column < max_column; ++column)
            {
                m_grid[row][column].set_randomResistance(mean, stddev);
            }
        }
    }

    //Recovery probability of the population
    void setRecoveryForAll(float const mean, float const stddev) //sets randomly the recovery probability of each Person
    {
        const int max_column = get_Length();
        const int max_row = get_Heigth();
        for(int row = 0; row < max_row; ++row)
        {
            for(int column = 0; column < max_column; ++column)
            {
                m_grid[row][column].set_randomRecoveryProbability(mean, stddev);
            }
        }
    }
    void setRecovery(float const prob) //sets the same recovery probability for every Person
    {
        const int max_column = get_Length();
        const int max_row = get_Heigth();

        for(int row = 0; row < max_row; ++row)
        {
            for(int column = 0; column < max_column; ++column)
            {
                m_grid[row][column].set_recoveryProbability(prob);
            }
        }
    }

    //Infection
    void set_InfectionProbability(float const prob) { m_infectionProbability = prob; }
    auto get_InfectionProbability() const { return m_infectionProbability; }

    void set_InfectionRadius(int const r) { m_infectionRadius = r; }
    auto get_InfectionRadius() const { return m_infectionRadius; }

    void setInfected(int const row, int const column) //sets the cell at coordinates [r][c] as Infected
    {
        m_grid[row][column].set_state(State::Infected);
    }
    void setInfected(Position const& cellPos) //sets the cell with the passed position as Infected
    {
        m_grid[cellPos.R][cellPos.C].set_state(State::Infected);
    }
    void setSomeInfected(int const num) //introduces some infected into the grid
    {
        int i;
        int j;
        const int max_columnIndex = get_Length() - 1;
        const int max_rowIndex = get_Heigth() - 1;
        for(int h = 0; h < num; ++h) //selects uniformly a row index and a column index
        {
            i = uniformIntValue(0, max_rowIndex);
            j = uniformIntValue(0, max_columnIndex);
            setInfected(i,j);
        }
    }

};

//returns a list of coordinates of all the cells which create a circumference of a certain radius centered in the selected cell
std::vector<Position> approxCircumference(Position const& centerCellPos, int const radius)
{
    int x;
    int y;
    int i; //related to rows
    int j; //related to columns
    int const width = 2 * radius + 1;
    int const circularSectors = 2 * 4 * width; // 4*width is the perimeter of the circumscribed square
    double const theta = (2 * pi) / circularSectors;  //a circle (range) is divided into a '8*width' number of equally wide angles (or sectors), going from 0 to 360°
    std::vector<Position> positions;
    
    for(int k = 0; k <= circularSectors; ++k)
    {
        x = radius * cos(k * theta); //casting to ints
        y = radius * sin(k * theta);

        i = centerCellPos.R + y;
        j = centerCellPos.C + x;

        Position newPosition{i,j};

        if(k > 0)
        {
            if(positions.back() == newPosition) { continue; } //ignores consequent duplicates
        }
        positions.push_back(newPosition); //executes at least once
    
    }
    return positions;
}

void createRange(std::vector< std::vector<Person> > & grid, Position const& centerCellPos, int const radius)
{
    int i;
    int j;
    int const max_row = grid.size();
    int const max_column = grid[0].size();
    
    std::vector<Position> cellsList = approxCircumference(centerCellPos, radius);

    for(int h = 0; h != cellsList.size(); ++h)
    {
        i = cellsList[h].R;
        j = cellsList[h].C;

        if(!(i < 0 || j < 0 || i >= max_row || j >= max_column))
        {
            grid[i][j].set_inRange(true); //if the selected cell is inside the grid, it gets set as in range
        }
        else
        {
            if(i < 0 || i >= max_row) //if the selected cell is (outside) above or below the grid, the cycle continues...
            {
                continue;
            }
            else //...if it is (outside) to its left or right, the closest cell in the same row becomes in range
            {
                if(j < 0) { grid[i][0].set_inRange(true); }
                if(j >= max_column) { grid[i][max_column - 1].set_inRange(true); }
            }
        } 
    }

    //NOTE: the approximated circle can be seen as a number of one-cell-tall horizontal rectangles stacked on top of eachother
    int leftEdge; //left-most square of a rectangle
    int rightEdge; //right-most square of a rectangle
    bool alreadyEncountered;
    
    //cycle that sets the inner cells within the circumference as inRange; this fills the entire circle
    for(i = centerCellPos.R - radius; i <= centerCellPos.R + radius; ++i)
    {
        alreadyEncountered = false;

        for(j = centerCellPos.C - radius; j <= centerCellPos.C + radius; ++j)
        {
            if(i < 0 || j < 0 || i >= max_row || j >= max_column || grid[i][j].get_inRange() == false)
            {
                continue;
            }
            else if(grid[i][j].get_inRange() == true) //extra check just to be sure -- could be omitted
            {  
                if(alreadyEncountered == false) //checks whether it's the first cell in range encountered
                {
                    leftEdge = j;
                    rightEdge = j;
                    alreadyEncountered = true;
                }
                else if(alreadyEncountered == true) //extra check just to be sure -- could be omitted
                {
                    if(j < leftEdge) { leftEdge = j; } //technically, this line is superflous since the loop is going from left to right regardless
                    if(j > rightEdge) { rightEdge = j; }
                }
            }
        }

        //if by cycling left to right through a row/rectangle a cell in range is found, it means that there is
        //a starting point from where to start 'filling'
        if(alreadyEncountered == true)
        {
            for(j = leftEdge; j <= rightEdge; ++j) //'fills' one horizontal rectangle at a time
            {
                if(grid[i][j].get_inRange() == false)
                {
                    grid[i][j].set_inRange(true);
                }
            }
        } 
    } 
}

//attempts to infect a cell by comparing the infection probability (corrected with the cell's resistance value) and
// a uniformly generated value
void tryInfect(Person & person, float const infectionProbability)
{
    auto const correctedInfProb = infectionProbability * (1 - person.get_resistance());
    auto const extractedResult = uniformRealValue(0., 1.);
    if(extractedResult <= correctedInfProb) //checks if the extracted value is part of the range whose upper bound is set by the final probability
    {
        person.set_state(State::Infected);
    }
}

//tries to infect Susceptible people around the selected Infected cell
//resets the inRange state of the cells it checks
void infectionAttempts(std::vector< std::vector<Person> > & grid, Position const& centerCell, int const radius, float const infectionProbability)
{
    int i; //rows
    int j; //columns
    int const max_row = grid.size();
    int const max_column = grid[0].size();

    for(i = centerCell.R - radius; i <= centerCell.R + radius; ++i)
    {
        for(j = centerCell.C - radius; j <= centerCell.C + radius; ++j)
        {
            if(i < 0 || j < 0 || i >= max_row || j >= max_column || grid[i][j].get_inRange() == false)
            {
                continue; //it's pointless to try infecting cells either outside of the grid, or out of range
            }
            else if(grid[i][j].get_inRange() == true) //could be omitted
            {
                if(grid[i][j].get_state() == State::Susceptible)
                {
                    tryInfect(grid[i][j], infectionProbability);
                }
                grid[i][j].set_inRange(false); //resetting inRange status
            }
        }
    }
}

void tryRecovery(Person & person)
{
    auto const extractedResult = uniformRealValue(0., 1.);
    if(extractedResult <= person.get_recoveryProbability()) //checks if the extracted value is part of the range whose upper bound is set by the final probability
    {
        person.set_state(State::Recovered);
    }
}

//makes the cell at coordinates cellPos switch place with another one situated either above, below, to the right, to the left of itself
void moveCell(std::vector< std::vector<Person> > & grid, Position const& cellPos)
{
    int const max_row = grid.size();
    int const max_column = grid[0].size();

    //Person const& movingPerson = grid[cellPos.R][cellPos.C];

    std::vector<Position> newAvailablePositions;

    //if the space above the cell is inside the grid, the 'north' direction is available
    if(!(cellPos.R - 1 < 0)) { newAvailablePositions.push_back({cellPos.R - 1, cellPos.C}); }

    //if the space to the right of  the cell is inside the grid, the 'east' direction is available
    if(!(cellPos.C + 1 >= max_column)) { newAvailablePositions.push_back({cellPos.R, cellPos.C + 1}); }

    //if the space below the cell is inside the grid, the 'south' direction is available
    if(!(cellPos.R + 1 >= max_row)) { newAvailablePositions.push_back({cellPos.R + 1, cellPos.C}); }

    //if the space to the left of the cell is inside the grid, the 'west' direction is available
    if(!(cellPos.C - 1 < 0)) { newAvailablePositions.push_back({cellPos.R, cellPos.C - 1}); }

    Position const newPosition = newAvailablePositions[uniformIntValue(0, newAvailablePositions.size() - 1)]; //chooses one of the available positions
    newAvailablePositions.clear(); //the vector served its purpose, thus it can be emptied -- I suppose this line can be omitted
    /*
    Person temp = grid[newPosition.R][newPosition.C]; //temporarily storing the Person which will be overwritten by movingPerson
    grid[newPosition.R][newPosition.C] = movingPerson; //'moving' movingPerson (our selected cell) to the new coordinates
    grid[cellPos.R][cellPos.C] = temp; //putting temp in movingPerson's original place
    */
    std::swap(grid[newPosition.R][newPosition.C], grid[cellPos.R][cellPos.C]);
}

void evolve(Board & board)
{
    int const max_column = board.get_Length(); //total height of the grid
    int const max_row = board.get_Heigth(); //total length of the grid

    auto someGrid = board.get_Grid(); //should create an instance initialised to the board's m_grid

    Position selectedCellPos;
    auto const infectionProbability = board.get_InfectionProbability();
    auto const radius = board.get_InfectionRadius();

    //indexes
    int i; //row
    int j; //column

    //movement
    for(i = 0; i < max_row; ++i) //cycles the rows of the grid
    {
        for(j = 0; j < max_column;) //cycles the columns of the grid
        {
            selectedCellPos = {i,j}; //sets the coordinates of the currently selected cell
            Person & person = someGrid[i][j];

            if(person.get_hasMoved() == false)
            {
                person.set_hasMoved(true);
                moveCell(someGrid, selectedCellPos);
            }
            else
            {
                ++j;
            }
        }
    }

    //recovery and infection
    for(i = 0; i < max_row; ++i)
    {
        for(j = 0; j < max_column; ++j)
        {
            selectedCellPos = {i,j};
            Person & person = someGrid[i][j];

            person.set_hasMoved(false); //resetting the hasMoved status

            //recovery attempt of the currently selected cell
            if(person.get_state() == State::Infected)
            {
                tryRecovery(person);
            }

            //attempts at infecting the cells around the one currently selected
            if(person.get_state() == State::Infected)
            {
                createRange(someGrid, selectedCellPos, radius); //sets inRange status
                infectionAttempts(someGrid, selectedCellPos, radius, infectionProbability); //resets inRange status
            }

        }
    }

    board.set_Grid(someGrid); //sets the modified grid into the board
}

void printHeadingToFile(Board const& board)
{
    std::time_t const& now = std::time(0);
    std::string fullDT = ctime(&now);

    std::ofstream output;
    output.open("SIR_OutputData.txt");
    output << fullDT << '\n';
    output << "Population number =\t\t" << board.get_Heigth() * board.get_Length() << '\n';
    output << "Infection probability =\t" << board.get_InfectionProbability() << '\n';
    output << "Infection radius =\t\t" << board.get_InfectionRadius() << '\n';
    output << "-----------------------------------------------------------------------------------" << '\n';
    output << "\tS, I, R" << '\n';
    output.close();
}

void appendDataToFile(Board const& board, int const iterationNumber)
{
    int const max_row = board.get_Heigth();
    int const max_column = board.get_Length();
    auto const& grid = board.get_Grid();
    int numSusceptible = 0;
    int numInfected = 0;
    int numRecovered = 0;
    for(int i = 0; i < max_row; ++i)
    {
        for(int j = 0; j < max_column; ++j)
        {
            switch(grid[i][j].get_state())
            {
                case State::Susceptible : { ++numSusceptible; }
                    break;
                case State::Infected : { ++numInfected; }
                    break;
                case State::Recovered : { ++numRecovered; }
                    break;
                default:
                    break;
            }
        }
    }
    std::ofstream output;
    output.open("SIR_OutputData.txt", std::ios_base::app); //appends data to file
    output << iterationNumber << ")  ";
    output << numSusceptible << ", " << numInfected << ", " << numRecovered << '\n';
    output.close();
}

int main()
{
    //Text font -- The font needs to be loaded here for the error pop-up window to work
    sf::Font font;
    bool const isFontLoaded = font.loadFromFile("./Fonts/consola.ttf"); //[font].ttf must be included in the directory the program runs in

    if (isFontLoaded == false)
    {
        std::cout << "Warning: font couldn't be loaded. Text will not be displayed.\n"
                     "Make sure to have a 'Fonts' folder containing a 'consola.ttf' font file in the same directory as the program\n";
    }

    try
    {   
        //
        // VARIABLES
        //
        const int grid_Length = 200; //length of the (rectangular) grid in terms of number of cells/people
        const int grid_Heigth = 150; //heigth of the (rectangular) grid in terms of number of cells/people

        const float Resistance_mean = 0.20; //value of the mean of a gaussian distribution
        const float Resistance_variability = 0.20; //value of the standard deviation of a gaussian distribution
        
    //    float const Recovery_mean = 0.30;
    //    float const Recovery_variability = 0.20;
        float const Recovery_probability = 0.15; //same value for everyone

        float const Infection_probability = 0.90;
        int const Infection_radius = 2; //radius of the infection range

    //    Board b{grid_Length, grid_Heigth, Infection_probability, Infection_radius};
        Board b{grid_Length, grid_Heigth, Infection_probability, Infection_radius, Recovery_probability};

        b.setResistanceForAll(Resistance_mean, Resistance_variability); //every cell gets initialised with a random resistance value
    //    b.setRecoveryForAll(Recovery_mean, Recovery_variability); //every cell gets initialised with a random recovery probability
    //    b.setRecovery(Recovery_probability);
    //    b.set_InfectionRadius(Infection_radius);

        b.setSomeInfected(5); //sets some infected cells in the grid

        //
        // GRAPHICS
        //
        
        const int cellSize = 5; //dimension of the cell's side in px
        const int horizontalWinSize = grid_Length * cellSize;
        const int verticalWinSize = grid_Heigth * cellSize;

        //Window
        auto const &screen = sf::VideoMode::getDesktopMode();
        //std::cout << "\nscreen width = " << screen.width << "\nscreen height = " << screen.height << '\n';
        if(horizontalWinSize >= screen.width - 15)
        {
            throw std::runtime_error{"Window width is too big: either decrease cellSize or grid_Length"};
        }
        if(verticalWinSize >= screen.height - 15)
        {
            throw std::runtime_error{"Window height is too big: either decrease cellSize or grid_height"};
        }
        if (cellSize < 1)
        {
            throw std::runtime_error{"cellSize cannot be less than 1"};
        }

        sf::RenderWindow window(sf::VideoMode(horizontalWinSize, verticalWinSize), "SIR Model simulation", sf::Style::Close);
        window.setPosition(sf::Vector2i(500, 100)); //change the position of the window (relatively to the desktop)
        window.setVerticalSyncEnabled(true);        //call it once, after creating the window
        window.setKeyRepeatEnabled(false);          //disables repeated "keyPressed" events

        // DRAWABLE OBJECTS //

        //Cell
        sf::RectangleShape cell{sf::Vector2f(cellSize, cellSize)}; //creates a square
        cell.setFillColor(sf::Color::White);

        // TEXT-RELATED //

        //Font is loaded at the start of main()
        
        //Caption background
        sf::RectangleShape captionBg{sf::Vector2f(230,55)};
        captionBg.setFillColor(sf::Color(130,130,130, 180));
        captionBg.setPosition(0,0);
        //auto bg_pos = captionBg.getPosition(); //needed to set the text's position relative to it
        

        //Caption text
        std::string controls_text{"Infect cell:\t  [LeftClick]\nToggle evolution: [RightClick]\nEvolution = "};
        int const captionOffset_x = 10;
        int const captionOffset_y = 5;

        sf::Text caption;
        caption.setFont(font);
        caption.setString(controls_text); //sets the text to display
        caption.setCharacterSize(12);
        caption.setFillColor(sf::Color::White);
        caption.setOutlineColor(sf::Color());
        //text.setStyle(sf::Text::Bold | sf::Text::Underlined); //if needed
        caption.setPosition(captionOffset_x, captionOffset_y); //its position is set relatively to the top-left corner of the window

        auto const &captionRect = caption.getLocalBounds(); //needed to check if the text should be displayed

        sf::Text OnOff;
        OnOff.setFont(font);
        //OnOff.setString(...) is down at the end
        OnOff.setCharacterSize(12);
        //setFillColor() is also down at the end
        OnOff.setOutlineColor(sf::Color());
        OnOff.setStyle(sf::Text::Bold);
        //OnOff.setPosition(bg_pos.x + captionOffset_x, bg_pos.y + captionOffset_y);
        OnOff.setPosition(captionOffset_x, captionOffset_y); //its position must be the same as the caption's

        //Watermark
        std::string watermark_text{"Daniel Michelin, July 2020"};

        sf::Text watermark;
        watermark.setFont(font);
        watermark.setString(watermark_text);
        watermark.setCharacterSize(12);
        watermark.setFillColor(sf::Color::White);
        watermark.setOutlineColor(sf::Color());

        auto const &watermarkRect = watermark.getLocalBounds(); //gets the parameters of the rectangle which corresponds to the text
        //apparently, the height of this rectangle is always less than the character size, which is dumb

        //std::cout << "\nwidth = " << watermarkRect.width << "\nheight = " << watermarkRect.height << '\n'; //for testing purposes

        int const compensatingHeight = 4;                                       // 2 is to have the text not cut off; the addittional 2 is to have commas not cut off
        watermark.setOrigin(watermarkRect.width - 1, watermarkRect.height - 1); //sets the bottom-right corner as the origin
        watermark.setPosition(horizontalWinSize - 1, verticalWinSize - 1 - compensatingHeight);   //anchors it to the bottom-right corner of the window

        // GAME LOOP //

        bool evolutionIsActive = false;
        int const evolutionPeriod = 300; //period of evolution in milliseconds
        sf::Clock clock;                 //the clock starts counting the time
        int iterationsCounter = 0;

        printHeadingToFile(b); //must be executed once
        
        while (window.isOpen()) //run the program as long as the window is open
        {
            sf::Time const &elapsedTime = clock.getElapsedTime();
            if (elapsedTime.asMilliseconds() >= evolutionPeriod)
            {
                if (evolutionIsActive == true)
                {
                    if(iterationsCounter == 0)
                    {
                        appendDataToFile(b, iterationsCounter); //outputs initial state
                    }
                    evolve(b);
                    ++iterationsCounter;
                    appendDataToFile(b, iterationsCounter);
                }
                clock.restart();
                //this way, elapsedTime will never increase indefinitely if the program were left to idle
            }

            // EVENT HANDLING //

            //check all the window's events that were triggered since the last iteration of the loop
            sf::Event event;
            while (window.pollEvent(event))
            {
                //all these switches are here for an easier addition of button behaviour later

                switch (event.type)
                {
                case sf::Event::Closed:
                {
                    window.close(); //closes the window
                }
                break;

                case sf::Event::MouseButtonPressed:
                {
                    switch (event.mouseButton.button)
                    {
                    case sf::Mouse::Right:
                    {
                        evolutionIsActive = !evolutionIsActive;
                        clock.restart(); //restarts the time
                    }
                    break;

                    case sf::Mouse::Left:
                    {
                        sf::Vector2i mousePointerCoords = sf::Mouse::getPosition(window); //get the local mouse position (relative to a window)

                        const int cellPos_x = floor(mousePointerCoords.x / cellSize); //column index of the cell
                        const int cellPos_y = floor(mousePointerCoords.y / cellSize); //row index of the cell

                        auto const& person = b.get_Grid()[cellPos_y][cellPos_x];
                        if(person.get_state() == State::Susceptible)
                        {
                            b.setInfected({cellPos_y, cellPos_x});
                        }
                    }
                    break;

                    default:
                        break;
                    }
                }
                break;

                /*
                case sf::Event::KeyPressed:
                {
                    switch (event.key.code)
                    {
                    case sf::Keyboard::Space
                    {
                        
                    }
                    break;

                    default:
                        break;
                    }
                }
                break;
                */

                default:
                    break;
                }
            }

            // DISPLAYING THE OBJECTS //

            window.clear(sf::Color::Black); //clears the window with black color. Default: window.clear();

            //cell coloring and display cycle
            //indexes: i--> row, j--> column
            auto const& grid = b.get_Grid();
            for (int i = 0; i < grid_Heigth; ++i) //cycles the rows of the grid
            {
                for (int j = 0; j < grid_Length; ++j) //cycles the columns of the grid
                {
                    auto cell_to_display = cell; //creates a white square
                    cell_to_display.setPosition(j * cellSize, i * cellSize); //sets its position in relation to the matrix indexes
                    switch(grid[i][j].get_state()) //changes its color according to its State
                    {
                        case State::Susceptible : { cell_to_display.setFillColor(sf::Color(0,200,0)); } //dark green
                            break;
                        case State::Infected : { cell_to_display.setFillColor(sf::Color(200,0,0)); } //dark red
                            break;
                        case State::Recovered : { cell_to_display.setFillColor(sf::Color(0,0,200)); } //dark blue
                            break;
                        default:
                            break;
                    }
                    window.draw(cell_to_display); //displays it (lights up corresponding pixels)
                }
            }

            //Draws all the text with its background if the following conditions are met
            //Note: given that the font gets loaded, all of these drawables are instantiated -- they simply may not get displayed

            if (isFontLoaded == true &&
                (captionRect.width + captionOffset_x <= horizontalWinSize) &&
                (captionRect.height + captionOffset_y <= verticalWinSize) &&
                (watermarkRect.width <= horizontalWinSize) &&
                (watermarkRect.height + compensatingHeight <= verticalWinSize))
            {
                //window.draw(watermarkBg);
                window.draw(watermark);
                window.draw(captionBg);
                window.draw(caption);

                if (evolutionIsActive == true)
                {
                    OnOff.setString("\n\n\t\t   ON");
                    OnOff.setFillColor(sf::Color::Green);
                }
                else
                {
                    OnOff.setString("\n\n\t\t   OFF");
                    OnOff.setFillColor(sf::Color::Red);
                }

                window.draw(OnOff);
            }

            window.display(); //end the current frame
        }
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
    }
    
}