module PauliTools

function hello(s)
    "Hello $s"
end

function bsp(a::AbstractArray{Bool}, b::AbstractArray{Bool})
    # circshift b by half its 1st dimension to emulate symplectic product
    mod.(a * circshift(b, size(b, 1)/2), 2)  # mod elements to base 2
end

end