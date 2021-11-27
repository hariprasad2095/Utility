function [handover_trig,measure]=meas_mode_selection(RSRP_serving,RSRP_nbr,thresh1,char_flag)
    measure=0;
    handover_trig=0;
    if char_flag == 0
        %% normal mode of measurment
        measure=1;
    else
    
        %% m_mode of measurement. Implement the new algorithm. 
        % if the measured RSRP is less than threshold, start measurement.
        % if the measured RSRP is greater than threshold, no measurement is
        % needed. 
        if thresh1 > RSRP_serving
            measure=1;
        end
    end

    %% Handover check
    % implementing second phase of threshold. 
    if RSRP_serving-RSRP_nbr > 3
        handover_trig=1;
    else
        handover_trig=0;
    end
end




