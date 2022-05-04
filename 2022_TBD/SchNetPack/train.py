import numpy as np
import os, sys
import schnetpack as spk
from schnetpack.datasets import QM9
import torch
from schnetpack.environment import TorchEnvironmentProvider

rcut = sys.argv[1]
feats = sys.argv[2]
gauss = sys.argv[3]
ints = sys.argv[4]
trade = sys.argv[5]
frac = int(sys.argv[6])

# Load the training data
data_dr = '.'
data = 'data.db'
db = spk.AtomsData(os.path.join(data_dr, data),
        available_properties=['energy', 'forces'],
        environment_provider=TorchEnvironmentProvider(int(rcut), torch.device('cuda')))

train, val, test = spk.train_test_split(
        data=db,
        num_train=int(0.80*len(db)),
        num_val=int(0.20*len(db)),
        split_file="split.npz",
    )


train_loader = spk.AtomsLoader(train, batch_size=32, shuffle=True)
val_loader = spk.AtomsLoader(val, batch_size=32)

means, stddevs = train_loader.get_statistics('energy', divide_by_atoms=True)
print('Mean energy / atom:', means['energy'])
print('Std. dev. energy / atom:', stddevs['energy'])

n_features = int(feats)

schnet = spk.representation.SchNet(
    n_atom_basis=n_features,
    n_filters=n_features,
    n_gaussians=int(gauss),
    n_interactions=int(ints),
    cutoff=int(rcut),
    cutoff_network=spk.nn.cutoff.CosineCutoff
)

output = spk.atomistic.Atomwise(
    n_in=n_features,
    property='energy',
    mean=means['energy'],
    stddev=stddevs['energy'],
    derivative='forces',
    negative_dr=True)

model = spk.AtomisticModel(representation=schnet, output_modules=output)

from torch.optim import Adam

# tradeoff
rho_tradeoff = float(trade)

# loss function
def loss(batch, result):
    diff_e = batch['energy']-result['energy']
    err_sq_e = torch.mean(diff_e**2)

    diff_f = batch['forces']-result['forces']
    err_sq_f = torch.mean(diff_f**2)

    return rho_tradeoff*err_sq_e + (1-rho_tradeoff)*err_sq_f

# build optimizer
optimizer = Adam(model.parameters(), lr=1e-4)

import schnetpack.train as trn
metrics = [spk.metrics.MeanAbsoluteError('energy'),
           spk.metrics.MeanAbsoluteError('forces')]

hooks = [
    trn.CSVHook(log_path='.', metrics=metrics),
    trn.ReduceLROnPlateauHook(
        optimizer,
        patience=20,
        factor=0.75,
        min_lr=1e-6,
        stop_after_min=True
    )
]

trainer = trn.Trainer(
    model_path='.',
    model=model,
    hooks=hooks,
    loss_fn=loss,
    optimizer=optimizer,
    train_loader=train_loader,
    validation_loader=val_loader,)

device = "cuda"
n_epochs = 200
trainer.train(device=device, n_epochs=n_epochs)

