for i in {1,4}; do
    for j in {13,44,101,105,388,1643,1932,2358,2394,2754}; do
        screen -dmS "session__{$i}_{$j}" sh -c "python3 simulation_inference.py $j $i; exec bash"
    done
done