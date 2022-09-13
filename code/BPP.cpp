#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <random>
#include <cmath>

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

protected:
    std::vector<Item*> m_items;
    int m_bin_capacity;

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
         brute_force(0);
    };

    void print_solution() override
    {
        int i = 1;
        for (const auto bin : m_best_bins) {
            std::cout << "Bin #" << i << ":" << std::endl
                      << *bin << std::endl;
            i++;
        }

        std::cout << "MinBinPacking = " << m_best_solution << std::endl;
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
    void crossover(Chromosome* other) 
    {
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

        std::cout << first_break_point_1 << " " << second_break_point_1 << std::endl;
        std::cout << first_break_point_2 << " " << second_break_point_2 << std::endl;

        m_chromosome.insert(m_chromosome.begin() + second_break_point_1,
                            cloned_bins_1.begin(), cloned_bins_1.end());

        other->m_chromosome.insert(other->m_chromosome.begin() + first_break_point_2,
                                   cloned_bins_2.begin(), cloned_bins_2.end());

        std::cout << "Insert clones: " << std::endl;
        std::cout << *this << std::endl;
        std::cout << *other << std::endl;
 
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
        std::cout << "Remove duplicates: " << std::endl;
        std::cout << *this << std::endl;
        std::cout << *other << std::endl;

        sort_items(removed_items_1);
        sort_items(removed_items_2);

        std::cout << "Removed items: " << std::endl;
        for (auto item : removed_items_1) {
            std::cout << item->weight() << " ";
        }
        std::cout << std::endl;
        for (auto item : removed_items_2) {
            std::cout << item->weight() << " ";
        }
        std::cout << std::endl;

        insert_removed_items(m_chromosome, removed_items_1);
        insert_removed_items(other->m_chromosome, removed_items_2);

        std::cout << "Insert removed: " << std::endl;
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
    GeneticAlgorithmBinPacking();
    virtual ~GeneticAlgorithmBinPacking();

    void solve() override
    {

    };

protected:

private: 
    
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


int main(int argc, const char *argv[])
{
    std::vector<Item*> items;
    int bin_capacity;

    load_data(argv[1], items, bin_capacity);

    std::cout << "Bin Capacity: " << bin_capacity << std::endl;
    std::cout << "items: [ ";
    for (const auto item : items) {
        std::cout << item->weight() << " ";
    }
    std::cout << "]" << std::endl;

    auto bin_packing = new BruteForceBinPacking(items, bin_capacity);

    bin_packing->solve();
    bin_packing->print_solution();

    auto chromosome = new Chromosome(items, bin_capacity);

    std::cout << *chromosome << std::endl;

    chromosome->mutation(1);
    
    std::cout << *chromosome << std::endl;

    for (auto item : items) {
        delete item;
    }

    delete chromosome;

    delete bin_packing;

    return 0;
}
