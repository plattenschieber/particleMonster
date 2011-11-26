class World_LC:World {
public:
    // data structures
    /// cells
    std::vector<Cell> cells;
    /// Number of cells in every dimension
    int cell_N[DIM];
    /// length of cells
    real cell_length[DIM];
    /// r_cut used for calculation of the cell length
    real cell_r_cut;
};
