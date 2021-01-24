using LightningArt
using Luxor
using LinearAlgebra
using ForwardDiff
@quickactivate "LightningArt"
import GeometryBasics

@svg begin
    circle(Luxor.Point(0, 0), 200, :stroke)
end

@svg begin
    line(Luxor.Point(0, 0), Luxor.Point(100,100), :stroke)
end


using LightGraphs, MetaGraphs


graph = DiGraph()
add_vertex!(graph)
collect(vertices(graph))

ggplot(graph) = gplot(graph, nodelabel=1:nv(graph), edgelabel=1:ne(graph))
add_edge!(graph, 1, 2)

mg = MetaDiGraph(graph, 1.0)

abstract type AbstractCurve end
struct BasicCurve end
struct FunctionCurve{F} <: AbstractCurve 
    curve::F
end
(curve::FunctionCurve)(t) = curve.curve(t)

add_vertex!(mg)
ggplot(mg)
set_props!(mg, Edge(1,2), Dict{Symbol, BasicCurve}(:curve=>BasicCurve()))
props(mg, Edge(1,2))
props(mg, Edge(1,3))
add_vertex!(mg)
add_edge!(mg, 1, 3)
ggplot(mg)

collect(edges(mg))
set_prop!(mg, 1, :position, Point2(0.0, -250.0))
set_prop!(mg, 2, :position, Point2(0.0, 10.0))
set_prop!(mg, 3, :position, Point2(10.0, 10.0))

function render_graph(metagraph, savename; render_vertices=true, svg=true)
    if svg
        @svg begin
            position_array = get_prop(metagraph, :position)
            if render_vertices
                for vert in vertices(metagraph)
                    circle(Luxor.Point(position_array[vert]...), 0.1, :stroke)
                end
            end
            for edge in edges(metagraph)
                parent_pos = position_array[edge.src]
                child_pos = position_array[edge.dst]
                line(Luxor.Point(child_pos...), Luxor.Point(parent_pos...), :stroke)
            end
        end
    else
        @png begin
            position_array = get_prop(metagraph, :position)
            if render_vertices
                for vert in vertices(metagraph)
                    circle(Luxor.Point(position_array[vert]...), 0.1, :stroke)
                end
            end
            for edge in edges(metagraph)
                parent_pos = position_array[edge.src]
                child_pos = position_array[edge.dst]
                line(Luxor.Point(child_pos...), Luxor.Point(parent_pos...), :stroke)
            end
        end
    end
end

render_graph(mg)
set_prop!(mg, :position, [GeometryBasics.Point2(0.0, 0.0)])

normsq(v) = dot(v, v)

function create_diffeq_graph(n, initial_array, initial_gradients, potential_function, α, β, γ; expansion_rate=0.05)
    ∇potential_function(p) = ForwardDiff.gradient(potential_function, p)
    metagraph = MetaDiGraph(n)
    position_array = initial_array
    potential_array = potential_function.(initial_array)
    ∇potential_array = initial_gradients
    ∇sign_array = [1.0, -1.0]
    set_prop!(metagraph, :position, position_array)
    topleft_corner = GeometryBasics.Point2(-300.0, -300.0)
    bottomright_corner = GeometryBasics.Point2(300.0, 300.0)
    difference_width = bottomright_corner .- topleft_corner
    for child_vertex in length(initial_array)+1:n
        # random parent, not very good
        #parent_vertex = rand(1:child_vertex-1)  # TODO: change this for more variety in sampling


        target_position = GeometryBasics.Point2(rand(2)) .* difference_width .+ topleft_corner
        target_potential = potential_function(target_position)


        #@show target_potential, potential_array[2]
        #@show target_potential .- potential_array[2]
        distance_function(j) = α .* norm(target_position .- position_array[j]) .+ β .* norm(target_potential .- potential_array[j])

        parent_vertex = 1
        if child_vertex == 3
            parent_vertex = 1
        elseif child_vertex == 4
            parent_vertex = 2
        else
            parent_distance, parent_vertex = findmin(map(distance_function, 1:child_vertex-1))
        end

        gradient = ∇potential_array[parent_vertex]
        #@show parent_vertex
        #@show gradient

        parent_position = position_array[parent_vertex]
        child_position = expansion_rate .* normalize(normalize(target_position .- parent_position) .+ γ .* ∇sign_array[parent_vertex] .* normalize(gradient)) .+ parent_position
        push!(position_array, child_position)
        push!(potential_array, potential_function(child_position))
        push!(∇potential_array, ∇potential_function(child_position))
        push!(∇sign_array, ∇sign_array[parent_vertex])


        add_edge!(metagraph, parent_vertex, child_vertex, :curve, FunctionCurve(t->t * (child_position .- parent_position) .* t .+ parent_position))

    end
    metagraph

end

function top_bottom_proto_rrt(n; expansion_rate=5.0)
    top_middle_point = GeometryBasics.Point2(0.0, -220.0)
    bottom_middle_point = GeometryBasics.Point2(0.0, 220.0)
    initial_array = [top_middle_point, bottom_middle_point]
    initial_gradients = [GeometryBasics.Point2(0.0, 1.0), GeometryBasics.Point2(0.0, -1.0)]
    create_diffeq_graph(n, initial_array, initial_gradients, identity, 1.0, 0.0; expansion_rate=expansion_rate)
end

function top_bottom_proto_rrt_potential(n; expansion_rate=5.0)
    charge_constant = 100.0
    potential_function(p) = charge_constant .* (-1.0 ./ (norm(p .- top_middle_point) + 1e-1) .+ 1.0 ./ (norm(p .- bottom_middle_point) + 1e-1))
    top_middle_point = GeometryBasics.Point2(0.0, -200.0)
    bottom_middle_point = GeometryBasics.Point2(0.0, 200.0)
    initial_array = [top_middle_point, bottom_middle_point]
    initial_gradients = [GeometryBasics.Point2(0.0, 1.0), GeometryBasics.Point2(0.0, 1.0)]
    create_diffeq_graph(n, initial_array, initial_gradients, potential_function, 0.3, 1.0, 1.2; expansion_rate=expansion_rate)
end

middle_point = GeometryBasics.Point2(0.0,245.0)

mg = top_bottom_proto_rrt(3;expansion_rate=1.5)
mg = top_bottom_proto_rrt_potential(60000;expansion_rate=1.0)
render_graph(mg; render_vertices=true, svg=false)