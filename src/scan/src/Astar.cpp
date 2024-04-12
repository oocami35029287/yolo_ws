#include "Astar.hpp"

// Grid2D operators
bool Astar::Grid2D::operator == (const Grid2D& grid_) {
    return (x == grid_.x && y == grid_.y);
}
Astar::Grid2D Astar::Grid2D::operator + (const Astar::Grid2D& grid_) {
    return {x + grid_.x, y + grid_.y};
}


// Node struct initialization
Astar::Node::Node(Grid2D grid_, Node *parent_) {
    parent = parent_;
    grid = grid_;
    g_val = h_val = 0;
    decision = -1;
}
ostream& operator<<(ostream& out, const Astar::Node& node) {
    out << "(x:" << node.grid.x << ", y:" << node.grid.y;
    return out;
}
int Astar::Node::get_score() {
    return g_val + h_val;
}


// Constructor
Astar::Solver::Solver() {
    directions_ = {
        {0, 1}, {1, 0}, {0, -1}, {1, 1},
        {1, -1},// {-1, 1}, {-1, -1}  {-1, 0},
    };
    // set_heuristic(&Heuristic::manhattan);
    set_heuristic(&Heuristic::euclidean);
    set_diagonal_move(true);
}

void Astar::Solver::set_diagonal_move(bool enable) {
    num_directions_ = (enable ? 5 : 4);
}

void Astar::Solver::set_heuristic(HeuristicFunc h_func) {
    h_func_ = std::bind(h_func, std::placeholders::_1, std::placeholders::_2);
}

Astar::Node* Astar::Solver::find_node(std::vector<Node*>& nodes_list, Grid2D grid) {
    for (auto node : nodes_list) {
        if (node->grid == grid) return node;
    }
    return nullptr;
}

bool Astar::Solver::solve_ros(nav_msgs::OccupancyGrid::ConstPtr map_msg_ptr, nav_msgs::Path::Ptr path, int start_idx, int goal_idx, double timeout_ms) {
    // map_ptr_ = nav_msgs::OccupancyGrid::Ptr(new nav_msgs::OccupancyGrid(*map_msg_ptr));
    map_ptr_ = map_msg_ptr;
    int map_width = (int)map_ptr_->info.width;
    int map_height = (int)map_ptr_->info.height;

    std::vector<Astar::Node*> open_set;
    std::vector<Astar::Node*> close_set;
    Astar::Node* cur_node = nullptr;

    // Convert (x, y) to occupancy grid index
    Grid2D start = {start_idx % map_width, start_idx / map_width};
    Grid2D goal = {goal_idx % map_width, goal_idx / map_width};
    open_set.push_back(new Node(start));

    // Set timeout
    ros::Duration timeout = ros::Duration(timeout_ms/1000);
    ros::Time begin_time = ros::Time::now();
    
    if(is_collision(goal)){
        path->header.stamp = ros::Time::now();
        // cout << "Invaild goal" << endl;
        return false;
    }

    flag_success_ = false;
    while(!open_set.empty()) {
        if(ros::Time::now() - begin_time > timeout){
            flag_success_ = false;
            break;
        }

        cur_node = *open_set.begin();
        for(auto node:open_set) {
            if(node->get_score() < cur_node->get_score())       // Find the lowest cost node
                cur_node = node;
        }

        if(cur_node->grid == goal) {                        // If get goal
            flag_success_ = true;
            break;
        }

        // Push current node to "visited nodes" list
        close_set.push_back(cur_node);
        open_set.erase(std::find(open_set.begin(), open_set.end(), cur_node));

        // Explore walkable node
        for(int i = 0; i < num_directions_; ++i) {  
            Grid2D tmp_grid(cur_node->grid + directions_[i]);
            if(find_node(close_set, tmp_grid) || is_collision(tmp_grid))
                continue;                                           // Skip visited node & skip wall 

            int total_cost =   cur_node->g_val + ((i<3)? 10 : 14);    // Balance cost between 4 & 8 directions
            
            Node* successor = find_node(open_set, tmp_grid);
            if(successor == nullptr){
                successor = new Node(tmp_grid, cur_node);           // Expand a new node from current node
                successor->g_val = total_cost;
                successor->h_val = h_func_(successor->grid, goal)+5*against_wall_cost(cur_node->grid);  // Calc the heuristic value
                successor->decision = i;
                open_set.push_back(successor);
            }else if(total_cost < successor->g_val){ 
                // Update non-visited but expanded node if find that has lower cost    
                successor->parent = cur_node;   
                successor->g_val = total_cost;
                successor->decision = i;
            }
        }
    } // while loop end

    if(flag_success_){
        while(cur_node != nullptr) {
            geometry_msgs::PoseStamped pose;
            pose.pose.position.x = cur_node->grid.x * map_ptr_->info.resolution \
                                    + map_ptr_->info.origin.position.x;
            pose.pose.position.y = cur_node->grid.y * map_ptr_->info.resolution \
                                    + map_ptr_->info.origin.position.y;
            path->poses.push_back(pose);
            cur_node = cur_node->parent;
        }
    }
    
    path->header.stamp = ros::Time::now();
    return flag_success_;
}

bool Astar::Solver::is_collision(Grid2D grid) {
    int width = map_ptr_->info.width;
    int height = map_ptr_->info.height;

    // Check the (x,y) is valid
    if(grid.x < 0 || grid.y < 0 || grid.x >= width || grid.y >= height)
        return true;

    // Check collision
    int map_idx = grid.y * width + (grid.x % width);

    // Due to (the second large gaussian value*100) || (prefer not to go to unknown space)
    if(map_ptr_->data[map_idx] > 80 || map_ptr_->data[map_idx] < 0)  
        return true;
    else
        return false;
}

int Astar::Solver::against_wall_cost(Grid2D grid) {
    int width = map_ptr_->info.width;
    int map_idx = grid.y * width + (grid.x % width);
    return map_ptr_->data[map_idx];
}

Astar::Grid2D Astar::Heuristic::getDelta(Grid2D source, Grid2D target) {
    return{ abs(source.x - target.x),  abs(source.y - target.y) };
}

int Astar::Heuristic::manhattan(Grid2D source, Grid2D target) {
    auto delta = std::move(getDelta(source, target));
    return static_cast<uint>(10 * (delta.x + delta.y));
}

int Astar::Heuristic::euclidean(Grid2D source, Grid2D target) {
    auto delta = std::move(getDelta(source, target));
    return static_cast<uint>(10 * sqrt(pow(delta.x, 2) + pow(delta.y, 2)));
}

int Astar::Heuristic::octagonal(Grid2D source, Grid2D target) {
    auto delta = std::move(getDelta(source, target));
    return 10 * (delta.x + delta.y) + (-6) * std::min(delta.x, delta.y);
}