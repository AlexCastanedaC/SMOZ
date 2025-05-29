module potentials

export softSphere

function softSphere(r, σ)
    return (σ / r)^6
end

# end module
end