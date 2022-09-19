#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include <cmath>
#include <chrono>

#define EXIT_CHECK(flag, msg) { \
    if (flag) { \
        std::cout << msg << std::endl; \
        std::exit(EXIT_FAILURE); \
    } \
}

//#define DEBUG_PRINT

class Item
{
public:
    Item(int weight) 
        : m_weight(weight)
    {

    };

    Item(const Item &other)
        : m_weight(other.m_weight)
    {

    }

    virtual ~Item() 
    {

    };

    int weight() { return m_weight; };

private:
    int m_weight;
};

class Bin
{
public:
    Bin(const int capacity) 
        : m_capacity(capacity)
        , m_weight(0)
        , m_items(0) 
    {

    };

    Bin(const Bin &other) 
        : m_capacity(other.m_capacity)
        , m_weight(other.m_weight)
        , m_items(other.m_items)
    {

    };

    virtual ~Bin() 
    {

    };

    bool put(Item *item) 
    {
        if (m_weight + item->weight() <= m_capacity) {
            m_items.push_back(item);
            m_weight += item->weight();
            return true;
        }

        return false;
    };

    void remove(Item *item) 
    {
        auto it = std::find(
                m_items.begin(), m_items.end(), 
                item
            );
        if (it == m_items.end()) {
            return;
        }
        m_weight -= (*it)->weight();
        m_items.erase(it);
    };

    void remove_all()
    {
        m_items.erase(m_items.begin(), m_items.end());
        m_weight = 0;
    }

    int number_of_items()
    {
        return m_items.size();
    };

    bool contains_some(const std::vector<Item*> &items)
    {
        for (auto item : m_items) {
            auto it = std::find(
                    items.begin(), items.end(), 
                    item
                );
            if (it != items.end()) {
                return true;
            }
        }

        return false;
    }

    int weight() { return m_weight; };
    std::vector<Item*> &items() { return m_items; };

protected:
    int m_capacity; 
    int m_weight;
    std::vector<Item*> m_items;

private:
    friend std::ostream &operator<<(std::ostream &os, const Bin &bin);
};

std::ostream &operator<<(std::ostream &os, const Bin &bin)
{
    os << "    Packing: [ ";
    for (auto item : bin.m_items) {
        os << item->weight() << " ";
    }
    os << "]" << std::endl
       << "    Weight: " << bin.m_weight
       << " (Capacity: " << bin.m_capacity << ")";
    return os;
}

class BinPacking
{
public:
    BinPacking(const std::vector<Item*> &items, const int bin_capacity)
        : m_items(items)
        , m_bin_capacity(bin_capacity)
        , m_duration()
        , m_solved(false)
    {

    }; 

    BinPacking(const BinPacking &other) 
        : m_items(other.m_items)
        , m_bin_capacity(other.m_bin_capacity)
    {

    };

    virtual ~BinPacking()
    {

    };

    virtual void solve() = 0;

    virtual void print_solution() = 0;

    virtual const int solution() = 0;

    virtual const std::string algorithm_type() = 0;

    constexpr int duration() const { return m_duration.count(); };

protected:
    std::vector<Item*> m_items;
    int m_bin_capacity;
    std::chrono::milliseconds m_duration;
    bool m_solved;
};

class BruteForceBinPacking: public BinPacking
{
public:
    BruteForceBinPacking(const std::vector<Item*> &items, 
                         const int bin_capacity)
        : BinPacking(items, bin_capacity)
    {
        for (auto item : items) {
            m_bins.push_back(new Bin(bin_capacity));
        }
        m_best_solution = m_bins.size();
    };

    virtual ~BruteForceBinPacking()
    {
        for (auto bin : m_bins) {
            delete bin;
        }
        for (auto bin : m_best_bins) {
            delete bin;
        }
    };

    void solve() override
    {
        auto start = std::chrono::high_resolution_clock::now();
        brute_force(0);
        auto end = std::chrono::high_resolution_clock::now();
        m_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        m_solved = true;
    };

    void print_solution() override
    {
        EXIT_CHECK(!m_solved, "Bin Packing is not solved!");

        int i = 1;
        for (const auto bin : m_best_bins) {
            std::cout << "Bin #" << i << ":" << std::endl
                      << *bin << std::endl;
            i++;
        }

        std::cout << "MinBinPacking = " << m_best_solution << std::endl;
    };

    const int solution() override
    {
        EXIT_CHECK(!m_solved, "Bin Packing is not solved!");

        return m_best_solution;
    };

    const std::string algorithm_type() override
    {
        return "BF";
    };

private:
    int count_filled_bins() 
    {
        int count = 0;
        for (const auto bin : m_bins) {
            if (bin->number_of_items() > 0) {
                count++;
            }
        }
        return count;
    };

    void update_best_bins(int new_solution)
    {
        if (m_best_solution > new_solution) {
            m_best_solution = new_solution;
            for (auto bin : m_best_bins) {
                delete bin;
            }
            m_best_bins.erase(m_best_bins.begin(), m_best_bins.end());
            for (auto bin : m_bins) {
                if (bin->number_of_items() > 0) {
                    m_best_bins.push_back(new Bin(*bin));
                }
            }
        }
    };

    void brute_force(int current_item_index) 
    {
        if (current_item_index >= m_items.size()) {
            update_best_bins(count_filled_bins());
            return;
        }

        const auto current_item = m_items[current_item_index];
        for (auto bin : m_bins) {
            if (bin->put(current_item)) {
                brute_force(current_item_index + 1);
                bin->remove(current_item);
            }
        }
    };

private:
    std::vector<Bin*> m_bins;
    std::vector<Bin*> m_best_bins;
    int m_best_solution;
};


class Chromosome
{
public:
    Chromosome(const std::vector<Item*>& items, const int bin_capacity) 
        : m_bin_capacity(bin_capacity)
        , m_chromosome()
        , m_fitness(0.0)
    {
        initialize(items);
    };

    virtual ~Chromosome()
    {

    };

private:
    void initialize(std::vector<Item*> items)
    {
        std::random_shuffle(items.begin(), items.end());

        Bin bin(m_bin_capacity);
        for (auto item : items) {
            if (!bin.put(item)) {
                m_chromosome.push_back(Bin(bin));
                bin.remove_all();
                bin.put(item);
            }
        }

        if (bin.number_of_items() > 0) {
            m_chromosome.push_back(Bin(bin));
        }
    };

public:
    void crossover(Chromosome* other, double probability) 
    {
        if (std::rand() / (double) RAND_MAX > probability) {
            return;
        }

        int n = m_chromosome.size();
        int m = other->m_chromosome.size();

        int first_break_point_1 = std::rand() % n;
        int second_break_point_1 = std::rand() % n;
        if (first_break_point_1 > second_break_point_1) {
            std::swap(first_break_point_1, second_break_point_1);
        }
        second_break_point_1++;

        int first_break_point_2 = std::rand() % m;
        int second_break_point_2 = std::rand() % m;
        if (first_break_point_2 > second_break_point_2) {
            std::swap(first_break_point_2, second_break_point_2);
        }
        second_break_point_2++;

        std::vector<Bin> cloned_bins_1;
        for (auto bin = other->m_chromosome.begin() + first_break_point_2;
             bin != other->m_chromosome.begin() + second_break_point_2;
             ++bin) {
            cloned_bins_1.push_back(Bin(*bin));
        };
        std::vector<Bin> cloned_bins_2;
        for (auto bin = m_chromosome.begin() + first_break_point_1;
             bin != m_chromosome.begin() + second_break_point_1;
             ++bin) {
            cloned_bins_2.push_back(Bin(*bin));
        };

        m_chromosome.insert(m_chromosome.begin() + second_break_point_1,
                            cloned_bins_1.begin(), cloned_bins_1.end());

        other->m_chromosome.insert(other->m_chromosome.begin() + first_break_point_2,
                                   cloned_bins_2.begin(), cloned_bins_2.end());

#ifdef DEBUG_PRINT
        std::cout << "Insert clones: " << std::endl;
        std::cout << *this << std::endl;
        std::cout << *other << std::endl;
#endif
 
        auto removed_items_1 = remove_duplicate_bins(
                cloned_bins_1, 
                m_chromosome, 
                second_break_point_1, 
                second_break_point_1 + (second_break_point_2 - first_break_point_2)
            );
        auto removed_items_2 = remove_duplicate_bins(
                cloned_bins_2, 
                other->m_chromosome, 
                first_break_point_2, 
                first_break_point_2 + (second_break_point_1 - first_break_point_1)
            );
#ifdef DEBUG_PRINT
        std::cout << "Remove duplicates: " << std::endl;
        std::cout << *this << std::endl;
        std::cout << *other << std::endl;
#endif

        sort_items(removed_items_1);
        sort_items(removed_items_2);

#ifdef DEBUG_PRINT
        std::cout << "Removed items: " << std::endl;
        for (auto item : removed_items_1) {
            std::cout << item->weight() << " ";
        }
        std::cout << std::endl;
        for (auto item : removed_items_2) {
            std::cout << item->weight() << " ";
        }
        std::cout << std::endl;
#endif

        insert_removed_items(m_chromosome, removed_items_1);
        insert_removed_items(other->m_chromosome, removed_items_2);

#ifdef DEBUG_PRINT
        std::cout << "Insert removed: " << std::endl;
        std::cout << *this << std::endl;
        std::cout << other << std::endl;
#endif
    };

    void mutation(double probability) 
    {
        int n = m_chromosome.size();

        std::vector<Item*> removed_items;
        auto bin = m_chromosome.begin();
        while (bin != m_chromosome.end()) {
            if (std::rand() / (double) RAND_MAX < probability) {
                for (auto item : bin->items()) {
                    removed_items.push_back(item);
                }
                bin = m_chromosome.erase(bin);
            } else {
                bin++;
            }
        }
        
        sort_items(removed_items);
        insert_removed_items(m_chromosome, removed_items);
    };

    double calculate_fitness(int k=2)
    {
        m_fitness = 0;
        for (auto bin : m_chromosome)
        {
            m_fitness += std::pow(bin.weight() / m_bin_capacity, k);
        }
        m_fitness /= m_chromosome.size();

        return m_fitness;
    };

    double fitness() { return m_fitness; };

    const std::vector<Bin> &chromosome() { return m_chromosome; };

private:
    void insert_removed_items(std::vector<Bin> &bins, 
                              const std::vector<Item*> &items) 
    {
        for (auto item : items) {
            bool flag = false;
            int i = 0;
            while (!flag) {
                if (i == bins.size()) {
                    Bin new_bin(m_bin_capacity);
                    new_bin.put(item);
                    bins.push_back(new_bin);
                    flag = true; 
                } else if (bins[i].put(item)) {
                    flag = true;
                }
                i++;
            }
        }
    };

    void sort_items(std::vector<Item*> &items) 
    {
        std::sort(
                items.begin(),
                items.end(),
                [] (auto item1, auto item2) {
                    return item1->weight() > item2->weight();
                }
        );
    };

    std::vector<Item*> remove_duplicate_bins(std::vector<Bin> &cloned_bins,
                                             std::vector<Bin> &bins, 
                                             int first_clone, 
                                             int last_clone)
    {
        std::vector<Item*> removed_items;

        std::vector<Item*> cloned_items;
        for (auto cloned_bin : cloned_bins) {
            for (auto cloned_item : cloned_bin.items()) {
                cloned_items.push_back(cloned_item);
            }
        }

        int i = 0;
        for (auto bin = bins.begin(); bin != bins.end(); bin++, i++) {
            if (i >= first_clone && i < last_clone) {
                continue;
            }

            if (bin->contains_some(cloned_items)) {
                for (auto item : bin->items()) {
                    auto it = std::find(
                            cloned_items.begin(), 
                            cloned_items.end(), 
                            item
                        );
                    if (it == cloned_items.end()) {
                        removed_items.push_back(item);
                    }
                }
                bin = bins.erase(bin);
                if (i < first_clone) {
                    first_clone--;
                    last_clone--;
                }
                bin--;
                i--;
                continue;
            }
        }

        return removed_items;
    };

private:
    std::vector<Bin> m_chromosome;
    double m_fitness;
    int m_bin_capacity;

private:
    friend std::ostream &operator<<(std::ostream &os, 
                                    Chromosome &chromosome);
};

std::ostream &operator<<(std::ostream &os, Chromosome &chromosome)
{
    std::cout << "Chromosome:" << std::endl;
    for (auto bin : chromosome.m_chromosome) {
        std::cout << bin << std::endl;
    }
    return os;
}

class GeneticAlgorithmBinPacking: public BinPacking
{
public:
    GeneticAlgorithmBinPacking(const std::vector<Item*> &items, 
                               const int bin_capacity,
                               const int population_size=100,
                               const int selection_size=50,
                               const int tournaments_size=5,
                               const int num_iterations=50,
                               const double crossover_probability=1.0,
                               const double mutation_probability=0.66)
        : BinPacking(items, bin_capacity)
        , m_population_size(population_size)
        , m_selection_size(selection_size)
        , m_tournament_size(tournaments_size)
        , m_num_iterations(num_iterations)
        , m_iter(0)
        , m_crossover_probability(crossover_probability)
        , m_mutation_probability(mutation_probability)
        , m_best_solution(nullptr)
    {

    };

    virtual ~GeneticAlgorithmBinPacking()
    {
        delete m_best_solution;
    };

    void solve() override
    {
        auto start = std::chrono::high_resolution_clock::now();
        genetic_alogorithm();
        auto end = std::chrono::high_resolution_clock::now();
        m_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        m_solved = true;
    };

    void print_solution() override
    {
        EXIT_CHECK(!m_solved, "Pin Packing is not solved!");

        const auto bins = m_best_solution->chromosome();

        int i = 1;
        for (const auto bin : bins) {
            std::cout << "Bin #" << i << ":" << std::endl
                      << bin << std::endl;
            i++;
        }

        std::cout << "MinBinPacking = " << bins.size() << std::endl;
    };

    const int solution() override
    {
        EXIT_CHECK(!m_solved, "Pin Packing is not solved!");

        return m_best_solution->chromosome().size();
    };

    const std::string algorithm_type() override
    {
        return "GA";
    };

private:

    Chromosome* selection(const std::vector<Chromosome*> &population) 
    {
        const int n = population.size();

        Chromosome* winner(nullptr);
        for (int i = 0; i < m_tournament_size; i++) {
            const int rand_index = rand() % n;
            if (winner == nullptr) {
                winner = population[rand_index];
            } else {
                if (winner->fitness() < population[rand_index]->fitness()) {
                    winner = population[rand_index];
                }
            }
        }

        return winner; 
    };

    void sort_chromosomes(std::vector<Chromosome*> &population) 
    {
        std::sort(
            population.begin(), population.end(), 
            [](auto chromosome1, auto chromosome2) {
                return chromosome1->fitness() > chromosome2->fitness();
            }
        );
    };

    bool stop_cond() 
    {
        if (m_iter < m_num_iterations) {
            return true;
        }

        return false;
    };

    void genetic_alogorithm()
    {
        std::vector<Chromosome*> population;
        for (int i = 0; i < m_population_size; i++) {
            Chromosome* chromosome = new Chromosome(m_items, m_bin_capacity);
            chromosome->calculate_fitness();
            population.push_back(chromosome);
        }
        sort_chromosomes(population);

        m_iter = 0;
        while(!stop_cond()) {
            for (int i = 0; i < m_selection_size; i++) {
                auto child1 = selection(population);
                auto child2 = selection(population);
                if (child1 == child2) {
                    continue;
                }
                child1->crossover(child2, m_crossover_probability);

                child1->mutation(m_mutation_probability);
                child2->mutation(m_mutation_probability);

                child1->calculate_fitness();
                child2->calculate_fitness();
            }
            sort_chromosomes(population);

            m_iter++;
        }
        
        m_best_solution = population[0];

        for (int i = 1; i < m_population_size; i++) {
            delete population[i];
        }
    };

private: 
    int m_population_size;
    int m_selection_size;
    int m_tournament_size;
    int m_num_iterations;
    int m_iter;

    double m_crossover_probability;
    double m_mutation_probability;

    Chromosome *m_best_solution;
};

void load_data(const std::string &file_name, 
               std::vector<Item*> &items, 
               int &bin_capacity)
{
    std::ifstream file;

    file.open(file_name, std::ifstream::in);

    int num_items;

    file >> num_items;
    file >> bin_capacity;

    for (int i = 0; i < num_items; i++) {
        int item;
        file >> item;
        items.push_back(new Item(item));
    }

    file.close();
}

void save_results(const std::string &file_name, 
                  const std::string &instance_name,
                  BinPacking* bin_packing)
{
    std::ofstream file;

    file.open(file_name, std::ofstream::app);

    file << instance_name << "," 
         << bin_packing->algorithm_type() << "," 
         << bin_packing->solution() << ","
         << bin_packing->duration() << std::endl;

    file.close();
}


int main(int argc, const char *argv[])
{
    std::string in_file(argv[1]);
    std::string out_file(argv[2]);

    auto instance_name = in_file.substr(in_file.find_last_of("/") + 1);
    instance_name.erase(instance_name.find(".BPP"));

    std::vector<Item*> items;
    int bin_capacity;

    load_data(in_file, items, bin_capacity);

#ifdef DEBUG_PRINT
    std::cout << "Bin Capacity: " << bin_capacity << std::endl;
    std::cout << "items: [ ";
    for (const auto item : items) {
        std::cout << item->weight() << " ";
    }
    std::cout << "]" << std::endl;
#endif

    if (argc > 3 && std::string(argv[3]) == std::string("-b")) {
        auto bin_packing = new BruteForceBinPacking(items, bin_capacity);

        bin_packing->solve();
#ifdef DEBUG_PRINT
        bin_packing->print_solution();
#endif

        save_results(out_file, instance_name, bin_packing);

        delete bin_packing;
    } else {
        auto bin_packing = new GeneticAlgorithmBinPacking(items, bin_capacity);

        bin_packing->solve();
#ifdef DEBUG_PRINT
        bin_packing->print_solution();
#endif

        save_results(out_file, instance_name, bin_packing);

        delete bin_packing;
    }

    for (auto item : items) {
        delete item;
    }

    return 0;
}
