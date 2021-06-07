module PauliTools

function hello(s)
    "Hello $s"
end

function bsp(a::AbstractArray{Bool}, b::AbstractArray{Bool})
    println("$(typeof(a)): $a")
    println("$(typeof(b)): $b")
    println("shifted b: ", circshift(b, size(b, 1)/2))
    mod.(a * circshift(b, size(b, 1)/2), 2)
end

end