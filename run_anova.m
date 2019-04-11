function [anova_out,model] = run_anova(tbl)
    model = fitlme(tbl,'y ~ 1 + groups*bins + (1|subject)','DummyVarCoding','reference');
    anova_out = anova(model,'DFMethod','satterthwaite');
end

