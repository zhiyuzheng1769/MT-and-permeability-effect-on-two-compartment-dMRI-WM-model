"""
repel_fixed_radius(radius, pos_in, repeatsize; maxiter=1000, repulsion_strength=0.01)
This function adds repulsion between randomly distributed cylinders generated by the base MCMRSimulator and obtains a more equally distributed configuration. It takes the fixed radius value, cylinder position output of random_positions_radii(), repetition block size as inputs, and outputs the refined cylinder positions. The repulsion_strength argument can be tuned to create configurations with differrent homogeneity.
""" 
function repel_fixed_radius(radius, pos_in, repeatsize; maxiter=1000, repulsion_strength=0.01)
    pos = deepcopy(pos_in)
    size = repeatsize
    
    function does_repulsed_have_a_collision(circle, index)
        if circle[1] > size/2 || circle[1] < -size/2
            return true
        elseif circle[2] > size/2 || circle[2] < -size/2
            return true
        end
        for i in 1:length(pos)
            if i != index
                other_circle = pos[i]
                a = 2*radius
                x = mod((circle[1] - other_circle[1] + size/2), size) - size/2
                y = mod((circle[2] - other_circle[2] + size/2), size) - size/2

                if a >= sqrt(x^2 + y^2)
                    return true
                end
            end
        end

        return false
    end

    # Update circle positions based on distance dependednt repulsion
    for a in 1:maxiter
        success_count=0
        for i in 1:length(pos)
            start_pos = pos[i]
            proposed_pos = pos[i]
            for j in 1:length(pos)
                if i != j
                    dx = mod((pos[j][1] - pos[i][1] + size/2), size) - size/2
                    dy =  mod((pos[j][2] - pos[i][2] + size/2), size) - size/2
                    dist = sqrt(dx^2 + dy^2)
                    force = (radius^2)*repulsion_strength / (dist)
                    proposed_pos = (proposed_pos[1] - force * dx/dist, proposed_pos[2] - force * dy/dist)
                end
            end
            # Abandon changes if there's overlap
            if does_repulsed_have_a_collision(proposed_pos, i)
                pos[i] = start_pos
            else
                pos[i] = proposed_pos
                success_count += 1
            end
        end
        # Stop if steady state/equilibrium is reached
        if success_count == 0
            println(a)
            break
        end
    end
        
    return pos
end

"""
repel_distributed_radius(radii, pos_in, repeatsize; maxiter=1000, repulsion_strength=0.01)
This function is the adaptation of the repel_fixed_radius() function to the case where the cylinder radii follows a Gamma distribution. The radius argument is now replaced with radii, a list of radius values corresponding to the elements in the pos_in argument recording positions of the cylinder centres.
""" 
function repel_distributed_radius(radii, pos_in, repeatsize; maxiter=1000, repulsion_strength=0.01)
    pos = deepcopy(pos_in)
    size = repeatsize
    
    function does_repulsed_have_a_collision(circle, index, pos)
        if circle[1] > size/2 || circle[1] < -size/2
            return true
        elseif circle[2] > size/2 || circle[2] < -size/2
            return true
        end
        for i in 1:length(pos)
            if i != index
                other_circle = pos[i]
                other_radius = radii[i]
                a = radii[index] + other_radius
                x = mod((circle[1] - other_circle[1] + size/2), size) - size/2
                y = mod((circle[2] - other_circle[2] + size/2), size) - size/2

                if a >= sqrt(x^2 + y^2)
                    return true
                end
            end
        end

        return false
    end
    

        # Update circle positions based on repulsion
    for a in 1:maxiter
        success_count=0
        for i in 1:length(pos)
            start_pos = pos[i]
            proposed_pos = pos[i]
            for j in 1:length(pos)
                if i != j
                    # ((abs(circle[1][1] - other_circle[1][1]) + size/2) % size) - size/2
                    dx = mod((pos[j][1] - pos[i][1] + size/2), size) - size/2
                    dy =  mod((pos[j][2] - pos[i][2] + size/2), size) - size/2
                    dist = sqrt(dx^2 + dy^2)
                    # force = (radius^2*radius^2)*repulsion_strength / (dist * mean_radius^2)
                    force = (radii[i]^2*radii[j]^2)*repulsion_strength / (dist*mean(radii)^2)
                    proposed_pos = (proposed_pos[1] - force * dx/dist, proposed_pos[2] - force * dy/dist)
                    # println(dx, dy)
                end
            end
            # Abandon changes if there's overlap
            if does_repulsed_have_a_collision(proposed_pos, i, pos)
                pos[i] = start_pos
            else
                pos[i] = proposed_pos
                success_count += 1
            end
        end
        if success_count == 0
            # println(a)
            break
        end
        # println(success_count)
    end
        
    return pos
end