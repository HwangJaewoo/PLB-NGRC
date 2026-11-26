%% Reproduces all results, the training times and the figures in order; one failed block does not stop the others (log in Run_all_log.txt)
function Run_all()
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, 'core'), genpath(fullfile(root, 'experiments')), fullfile(root, 'timing'), fullfile(root, 'figures', 'code'));
    blocks = {
        'Lorenz-96 system', @() Run_model_system('Lorenz96')
        'coupled Ikeda map', @() Run_model_system('Ikeda')
        'Thomas system', @() Run_model_system('Thomas')
        'Kuramoto-Sivashinsky equation', @KS_run
        'Z500 along 50N', @Z500_1D_run
        'Z500 over the latitude band', @Z500_2D_run
        'training time', @Measure_training_time
        'figures', @Make_figures
        };
    logf = fullfile(root, 'Run_all_log.txt');
    Note(logf, sprintf('==== start %s', datestr(now, 'yyyy-mm-dd HH:MM')));
    for i = 1:size(blocks, 1)
        Note(logf, sprintf('[%d/%d] %s: start %s', i, size(blocks, 1), blocks{i, 1}, datestr(now, 'HH:MM')));
        try
            blocks{i, 2}();
            Note(logf, sprintf('[%d/%d] %s: done %s', i, size(blocks, 1), blocks{i, 1}, datestr(now, 'HH:MM')));
        catch ME
            Note(logf, sprintf('[%d/%d] %s: FAILED %s: %s', i, size(blocks, 1), blocks{i, 1}, datestr(now, 'HH:MM'), ME.message));
        end
        diary off;
    end
    Note(logf, sprintf('==== end %s', datestr(now, 'yyyy-mm-dd HH:MM')));
end

function Make_figures()
    L = dir(fullfile(fileparts(mfilename('fullpath')), 'figures', 'code', 'Fig_*.m'));
    for i = 1:numel(L)
        Run_script(fullfile(L(i).folder, L(i).name));
        close all;
    end
end

function Run_script(f)
    run(f);
end

function Note(logf, s)
    fprintf('%s\n', s);
    fid = fopen(logf, 'a');
    if fid > 0
        fprintf(fid, '%s\n', s);
        fclose(fid);
    end
end
